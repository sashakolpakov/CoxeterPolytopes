"""Discrete source-recipe mapping for the frozen H5 proposal banks."""
FIXED=['seed2','compact1','compact2','compact3']
COMMON=['seed1']+FIXED
T=['seed6','crosspair0','nextH4','compact0'];BASES=['wall17','facePrism']
BANKS=[('fresh_h9_higher_face_20261013_scout','proposals'),('fresh_h9_higher_face_20261013_branch_scout','bank'),('fresh_h9_higher_face_20261013_profile_scout','candidates'),('fresh_h9_higher_face_20261013_finite_truncation_scout','candidates'),('fresh_h9_higher_face_20261013_truncation_scout','candidates'),('fresh_h9_h5_closure_20261013_conic_scout','candidates'),('fresh_h9_h5_closure_20261013_conic_scout_stage1','candidates'),('fresh_h9_h5_closure_20261013_conic_scout_stage2','candidates'),('fresh_h9_h5_closure_20261013_endpoint_scout','accepted'),('fresh_h9_h5_closure_next_20261013_scout','candidates'),('fresh_h9_h5_closure_next_20261013_stage1_scout','candidates'),('fresh_h9_h5_closure_next_20261013_stage2_scout','candidates')]
def edge(ci):return[n for i,n in enumerate(T)if i!=ci//2]+[BASES[ci%2]]+FIXED
def canonical(tag,j,row):
 if tag=='fresh_h9_higher_face_20261013_scout':return['crosspair0','nextH4','compact0','facePrism']+[n for n in COMMON if n!=row['free']]+['wall17'],[2]*8+[row['label']]
 if tag=='fresh_h9_higher_face_20261013_branch_scout':
  if j==0:return edge(1)+['wall17'],[2]*8+[16]
  if j<6:return edge(row['corner'])+['roof1'if j==4 else'roof0'],row.get('profile',[2]*8)+[row['anchor_label']]
  return row['edge']+['roof'+str(row['anchor_existing'])],row['profile']+[row['anchor_label']]
 if tag=='fresh_h9_higher_face_20261013_profile_scout':return row['edge']+['roof'+str(row['anchor_existing'])],row['profile']+[row['anchor_label']]
 if tag=='fresh_h9_higher_face_20261013_finite_truncation_scout':return row['right_fixed']+row['incident'],[2]*4+row['profile']
 if tag=='fresh_h9_higher_face_20261013_truncation_scout':return row['right']+[row['anchor']],[2]*8+[row['label']]
 ss=FIXED+row.get('four',row.get('edge'))+[row['anchor']];ll=[2]*4+row['profile']+[row['label']]
 if tag.endswith('_stage1_scout')and'_closure_next_'in tag:ss=[n.replace('endpoint','alternate_endpoint')if n.startswith('endpoint')else n for n in ss]
 return ss,ll
