"""Candidate reconstruction over Q or Q(sqrt(d)), checked at sqrt(d)>0.

LLL proposes kernels only. None is a universal constraint. Exact affine
projection and an exact PSD test decide whether a proposed point is accepted.
"""
from fractions import Fraction
from math import isqrt


class Field:
    def __init__(self,d=0):
        from flint import fmpq,ctx
        ctx.threads=1
        if d and (d<2 or isqrt(d)**2==d):raise ValueError('Use a positive nonsquare radicand')
        self.d=d;self.Q=fmpq;self.zero=(fmpq(0),fmpq(0));self.one=(fmpq(1),fmpq(0))

    def element(self,a=0,b=0):
        result=self.Q(str(a)),self.Q(str(b))
        if not self.d and result[1]:raise ValueError('A rational point cannot have a radical component')
        return result
    def add(self,x,y):return x[0]+y[0],x[1]+y[1]
    def neg(self,x):return -x[0],-x[1]
    def sub(self,x,y):return x[0]-y[0],x[1]-y[1]
    def mul(self,x,y):
        if not x[1] and not y[1]:return x[0]*y[0],self.Q(0)
        return x[0]*y[0]+self.d*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
    def scale(self,c,x):return c*x[0],c*x[1]
    def div(self,x,y):
        if not y[1]:return x[0]/y[0],x[1]/y[0]
        norm=y[0]*y[0]-self.d*y[1]*y[1]
        if not norm:raise ZeroDivisionError('Zero field denominator')
        return ((x[0]*y[0]-self.d*x[1]*y[1])/norm,
                (x[1]*y[0]-x[0]*y[1])/norm)
    def sign(self,x):
        a,b=x
        sa=(a>0)-(a<0);sb=(b>0)-(b<0)
        if not b:return sa
        if not self.d:raise ValueError('A nonrational coefficient was used over Q')
        if not a or sa==sb:return sb
        difference=a*a-self.d*b*b
        if not difference:raise ValueError('A nonsquare-field norm unexpectedly vanished')
        return sa if difference>0 else sb
    def floating(self,x):return float(x[0])+float(x[1])*(self.d**.5)
    def text(self,x):return [str(x[0]),str(x[1])]


def evaluate(row,point,field):
    a=row.get(-1,field.Q(0));b=field.Q(0)
    for col,c in row.items():
        if col>=0:
            a+=c*point[col][0]
            if point[col][1]:b+=c*point[col][1]
    return a,b


def psd(matrix,field):
    M=[list(row) for row in matrix];pivots=[]
    for i in range(len(M)):
        p=M[i][i];s=field.sign(p);pivots.append(p)
        if s<0:raise ValueError('Negative exact field LDL pivot at '+str(i))
        if not s:
            if any(M[i][j]!=field.zero for j in range(i+1,len(M))):
                raise ValueError('Nonzero exact row at zero field LDL pivot '+str(i))
            continue
        for j in range(i+1,len(M)):
            factor=field.div(M[j][i],p)
            for k in range(j,len(M)):
                M[j][k]=M[k][j]=field.sub(M[j][k],field.mul(factor,M[i][k]))
    return pivots


def check_point(forms,blocks,point,field,progress=None):
    checked=0;singular=0
    for block in blocks:
        n=block['size'];start=block['start']
        if n==1:
            value=evaluate(forms[start],point,field)
            if field.sign(value)<0:
                return dict(valid=False,block=block['name'],block_index=checked,
                    reason='Negative exact scalar',value=field.text(value),value_float=field.floating(value))
            singular+=value==field.zero
        else:
            M=[[field.zero for _ in range(n)] for _ in range(n)]
            pos=start
            for i in range(n):
                for j in range(i,n):
                    M[i][j]=M[j][i]=evaluate(forms[pos],point,field);pos+=1
            try:pivots=psd(M,field)
            except ValueError as error:
                return dict(valid=False,block=block['name'],block_index=checked,reason=str(error))
            singular+=any(p==field.zero for p in pivots)
        checked+=1
        if progress and checked%50000==0:progress(stage='exact_PSD_check',cones=checked)
    return dict(valid=True,cones_checked=checked,singular_cones=singular,
        exact_arithmetic='Q' if not field.d else 'Q(sqrt('+str(field.d)+')); sqrt(d)>0')


def rational_point(values,field,maximum_denominator):
    """Use one grid, avoiding unrelated denominators across a large SDP.

    This is candidate rounding only: it adds no information to binary64 data.
    Algebraic dependencies, when present, are imposed by exact projection.
    """
    point=[]
    for value in values:
        point.append((field.Q(round(float(value)*maximum_denominator),maximum_denominator),field.Q(0)))
    return point


def kernel_vectors(matrix,d=0,tolerance=1e-10,maximum_coefficient=128):
    """LLL proposals in actual matrix coordinates, never rounded eigenvectors."""
    import numpy as np
    from flint import fmpz_mat
    M=np.asarray(matrix,dtype=float);n=len(M)
    eigen,N=np.linalg.eigh(M)
    bound=tolerance*max(1.,float(np.max(np.abs(eigen))))
    nullity=int(np.count_nonzero(np.abs(eigen)<=bound))
    if not nullity:return []
    if nullity==n:
        return [[(int(i==j),0) for i in range(n)] for j in range(n)]
    # Direct integer combinations a+b*sqrt(d) are tested against the whole
    # matrix. This avoids mistaking a random eigenbasis for algebraic data.
    scale=10**10
    rows=n if not d else 2*n
    lattice=[]
    normalization=max(1.,float(np.max(np.abs(M))))
    for i in range(rows):
        row=[int(i==j) for j in range(rows)]
        multiplier=1. if i<n else d**.5
        row += [int(round(scale*multiplier*M[i%n,j]/normalization)) for j in range(n)]
        lattice.append(row)
    reduced=fmpz_mat(lattice).lll()
    vectors=[];numerical=[]
    for i in range(reduced.nrows()):
        integers=[int(reduced[i,j]) for j in range(rows)]
        height=max(map(abs,integers),default=0)
        if not height or height>maximum_coefficient:continue
        vector=[(integers[j],0 if not d else integers[n+j]) for j in range(n)]
        floating=np.array([a+b*(d**.5) for a,b in vector])
        norm=float(np.max(np.abs(floating)))
        if norm<1e-3*height:continue  # Reject tiny Diophantine approximations to zero.
        if float(np.max(np.abs(M@floating)))>100*bound*max(1.,norm):continue
        proposed=numerical+[floating]
        if np.linalg.matrix_rank(np.asarray(proposed),tol=1e-7)<=len(numerical):continue
        vectors.append(vector);numerical.append(floating)
        if len(vectors)==nullity:break
    return vectors


def add_kernel_equations(forms,block,vectors,field):
    n=block['size'];pos=block['start'];M=[[None for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i,n):M[i][j]=M[j][i]=forms[pos];pos+=1
    answer=[]
    for vector in vectors:
        for row in M:
            result={}
            for q,(a,b) in zip(row,vector,strict=True):
                if not a and not b:continue
                for col,c in q.items():
                    old=result.get(col,field.zero)
                    value=(old[0]+a*c,old[1]+b*c)
                    if value!=field.zero:result[col]=value
                    elif col in result:del result[col]
            if result:answer.append(result)
    return answer


def split_equations(equations,variables,field):
    """Restriction of scalars to two rational systems for a quadratic field."""
    answer=[]
    for row in equations:
        real={};imag={}
        for col,(a,b) in row.items():
            if col<0:
                if a:real[-1]=a
                if b:imag[-1]=b
            else:
                if a:real[col]=a
                if b:real[variables+col]=field.d*b
                if b:imag[col]=b
                if a and field.d:imag[variables+col]=a
        if real:answer.append(real)
        if imag:answer.append(imag)
    return answer


def project(equations,approximate,field,denominator=10**8,maximum_terms=12000000,maximum_bits=4096,progress=None):
    """Sparse exact affine solve with pointwise triangular backsubstitution."""
    from collections import Counter
    Q=field.Q;variables=len(approximate)
    rows=split_equations(equations,variables,field)
    frequency=Counter(k for row in rows for k in row if k>=0)
    basis={};order=[];terms=0
    for count,incoming in enumerate(sorted(rows,key=len)):
        if progress and count%10000==0:
            progress(stage='exact_candidate_elimination',rows=count,rank=len(basis),basis_terms=terms)
        row=dict(incoming)
        while True:
            old=next((k for k in row if k in basis),None)
            if old is None:break
            factor=row[old]
            for col,c in basis[old].items():
                row[col]=row.get(col,Q(0))-factor*c
                if not row[col]:del row[col]
            if len(row)>20000:raise ValueError('Candidate projection row-fill guard')
        names=[k for k in row if k>=0]
        if not names:
            if row:raise ValueError('Proposed candidate face is exactly inconsistent')
            continue
        pivot=min(names,key=lambda k:(frequency[k],-k))
        coefficient=row[pivot]
        normalized={k:v/coefficient for k,v in row.items()}
        if any(max(int(abs(v.numerator)).bit_length(),int(v.denominator).bit_length())>maximum_bits
               for v in normalized.values()):raise ValueError('Candidate projection coefficient-height guard')
        basis[pivot]=normalized;order.append(pivot);terms+=len(normalized)
        if terms>maximum_terms:raise ValueError('Candidate projection total-fill guard')
    initial=rational_point(approximate,field,denominator)
    values={k:v[0] for k,v in enumerate(initial)}
    if field.d:values.update({variables+k:Q(0) for k in range(variables)})
    for pivot in reversed(order):
        values[pivot]=-sum((c*(Q(1) if col<0 else values[col])
            for col,c in basis[pivot].items() if col!=pivot),Q(0))
        if max(int(abs(values[pivot].numerator)).bit_length(),int(values[pivot].denominator).bit_length())>maximum_bits:
            raise ValueError('Candidate field-point coefficient-height guard')
    if any(sum((c*(Q(1) if col<0 else values[col]) for col,c in row.items()),Q(0)) for row in rows):
        raise ValueError('A proposed kernel equation failed exact affine projection')
    point=[(values[k],values[variables+k] if field.d else Q(0)) for k in range(variables)]
    return point,dict(proposed_field_equations=len(equations),rational_equations=len(rows),
        exact_projection_rank=len(basis),triangular_basis_terms=terms,
        rational_free_denominator_bound=denominator)
