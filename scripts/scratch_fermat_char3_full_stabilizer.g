f := GF(9);
G := GU(4,3);
e := IdentityMat(4,f);
U := [e[1],e[2]];
lines := Orbit(G,U,OnSubspacesByCanonicalBasis);
Print("line_count=",Length(lines),"\n");

gram := List([1..Length(lines)], function(i)
  return List([1..Length(lines)], function(j)
    if i=j then
      return -2;
    elif RankMat(Concatenation(lines[i],lines[j]))=3 then
      return 1;
    else
      return 0;
    fi;
  end);
end);

LoadPackage("digraphs");
Read("/private/tmp/fermat_char3_archived_graph.g");
abstractAdjacency := List([1..112],i->Filtered([1..112],j->gram[i][j]=1));
Dabstract := Digraph(abstractAdjacency);
Darchived := Digraph(archivedAdjacency);
lineIso := IsomorphismDigraphs(Dabstract,Darchived);
if lineIso=fail then Error("abstract and archived line graphs are not isomorphic"); fi;
Print("line_graph_isomorphism=",ListPerm(lineIso,112),"\n");

N := gram{[1..22]}{[1..22]};
Print("first22_rank=",RankMat(N)," det=",DeterminantMat(N),"\n");

basis := [11,30,34,36,41,44,45,49,56,69,70,72,77,80,81,91,93,96,100,106,107,108];
N := gram{basis}{basis};
Print("basis=",basis," det=",DeterminantMat(N),"\n");

hom := ActionHomomorphism(G,lines,OnSubspacesByCanonicalBasis,"surjective");
P := Image(hom);
Print("PGU_order=",Size(P)," generators=",Length(GeneratorsOfGroup(P)),"\n");

FrobSubspace := function(B)
  local C;
  C := List(B,row->List(row,x->x^3));
  TriangulizeMat(C);
  return C;
end;

frob := PermList(List([1..Length(lines)],i->Position(lines,FrobSubspace(lines[i]))));
K := Group(Concatenation(GeneratorsOfGroup(P),[frob]));
Print("chamber_group_order=",Size(K)," frob_order=",Order(frob),"\n");

TransportToArchive := p->lineIso^-1*p*lineIso;
Karchive := Group(List(Concatenation(GeneratorsOfGroup(P),[frob]),TransportToArchive));
Print("archived_chamber_group_order=",Size(Karchive),"\n");

officialFrobImages := List([1..112], function(i)
  return Position(archivedLineClasses,archivedLineClasses[i]*archivedFrobenius);
end);
if fail in officialFrobImages then Error("FQAF does not permute the archived lines"); fi;
officialFrobPerm := PermList(officialFrobImages);
Print("official_frob_in_group=",officialFrobPerm in Karchive,
  " official_frob_order=",Order(officialFrobPerm),"\n");

ArchivedNSMatrix := function(p)
  local M;
  M := List(archivedNSBasis,i->archivedLineClasses[i^p]);
  if M*archivedGramS*TransposedMat(M)<>archivedGramS then
    Error("archived NS matrix does not preserve FQGramS");
  fi;
  return M;
end;

Ninv := N^-1;
NSMatrix := function(p)
  local M, i, img, ints, coords;
  M := [];
  for i in [1..22] do
    img := basis[i]^p;
    ints := gram[img]{basis};
    coords := ints*Ninv;
    if not ForAll(coords,IsInt) then
      Error("non-integral NS coordinates");
    fi;
    Add(M,coords);
  od;
  if M*N*TransposedMat(M)<>N then
    Error("NS matrix does not preserve Gram matrix");
  fi;
  return M;
end;

Print("generator_fixed_ranks=");
for p in Concatenation(GeneratorsOfGroup(P),[frob]) do
  M := NSMatrix(p);
  Print(22-RankMat(M-IdentityMat(22))," ");
od;
Print("\n");
Print("archived_generator_fixed_ranks=");
for p in GeneratorsOfGroup(Karchive) do
  M := ArchivedNSMatrix(p);
  Print(22-RankMat(M-IdentityMat(22))," ");
od;
Print("\n");

start := Runtime();
classes := ConjugacyClasses(K);
Print("class_count=",Length(classes)," class_runtime_ms=",Runtime()-start,"\n");
hist := [];
high := [];
highRecords := [];
involutionClasses := [];
highClassObjects := [];
highClassRanks := [];
highClassKeys := [];
for class_index in [1..Length(classes)] do
  c := classes[class_index];
  p := Representative(c);
  M := NSMatrix(p);
  fixedrank := 22-RankMat(M-IdentityMat(22));
  if Order(p)=2 then
    Add(involutionClasses,[Size(c),fixedrank,p in P,Number([1..112],i->i^p=i)]);
  fi;
  if not IsBound(hist[fixedrank+1]) then hist[fixedrank+1] := 0; fi;
  hist[fixedrank+1] := hist[fixedrank+1]+Size(c);
  if fixedrank>10 then
    Add(high,[class_index,Size(c),Order(p),fixedrank,CycleStructurePerm(p),p in P,
      Number([1..112],i->i^p=i)]);
    q := TransportToArchive(p);
    Add(highRecords,rec(
      class_index_in_this_run := class_index,
      class_size := Size(c),
      centralizer_order := Size(K)/Size(c),
      element_order := Order(p),
      fixed_rank := fixedrank,
      hyperbolic_dimension := fixedrank-1,
      in_projective_unitary_subgroup := p in P,
      fixed_archived_lines := Number([1..112],i->i^q=i),
      cycle_structure := CycleStructurePerm(q),
      representative_archived_line_images := ListPerm(q,112),
      representative_archived_ns_matrix := ArchivedNSMatrix(q)
    ));
    if Order(p)>1 then
      Add(highClassObjects,c);
      Add(highClassRanks,fixedrank);
      Add(highClassKeys,[p in P,Number([1..112],i->i^p=i)]);
    fi;
  fi;
od;
Print("fixed_rank_element_histogram=",hist,"\n");
Print("high_classes=",high,"\n");
Print("involution_classes=",involutionClasses,"\n");
Print("frob_class=",PositionProperty(classes,c->frob in c),"\n");
Print("official_frob_fixed_rank=",22-RankMat(archivedFrobenius-IdentityMat(22)),"\n");

HighFixedRank := function(x)
  local key, pos;
  if IsOne(x) then return 22; fi;
  if Order(x)<>2 then return fail; fi;
  key := [x in P,Number([1..112],i->i^x=i)];
  pos := Position(highClassKeys,key);
  if pos=fail then return fail; fi;
  return highClassRanks[pos];
end;

highElementBlocks := List(highClassObjects,AsList);
pairSummary := [];
eligibleDistinctPairs := [];
maxDistinctPairCommonFixedRank := 0;
for i in [1..Length(highClassObjects)] do
  a := Representative(highClassObjects[i]);
  fa := highClassRanks[i];
  localPairs := [];
  eligible := [];
  for j in [1..Length(highClassObjects)] do
    fb := highClassRanks[j];
    for b in highElementBlocks[j] do
      if b<>a then
        fp := HighFixedRank(a*b);
        if fp<>fail then
          commonrank := (fa+fb+fp-22)/2;
          Add(localPairs,[j,fp,commonrank]);
          if commonrank>10 then Add(eligible,b); fi;
          if commonrank>maxDistinctPairCommonFixedRank then
            maxDistinctPairCommonFixedRank := commonrank;
          fi;
        fi;
      fi;
    od;
  od;
  Add(pairSummary,[i,fa,Collected(localPairs)]);
  Add(eligibleDistinctPairs,eligible);
od;
Print("high_class_keys=",List([1..Length(highClassKeys)],i->[i,highClassRanks[i],
  Size(highClassObjects[i]),highClassKeys[i]]),"\n");
Print("distinct_high_pair_summary=",pairSummary,"\n");
Print("max_distinct_pair_common_fixed_rank=",maxDistinctPairCommonFixedRank,"\n");

GroupCommonFixedRank := function(H)
  local mats;
  mats := List(GeneratorsOfGroup(H),
    x->ArchivedNSMatrix(TransportToArchive(x))-IdentityMat(22));
  if Length(mats)=0 then return 22; fi;
  return 22-RankMat(Concatenation(mats));
end;

SubgroupProfile := function(H)
  return Collected(List(Elements(H),HighFixedRank));
end;

ConjugacyDeduplicateSubgroups := function(candidates)
  local reps, H;
  reps := [];
  for H in candidates do
    if ForAll(reps,J->RepresentativeAction(K,H,J)=fail) then
      Add(reps,H);
    fi;
  od;
  return reps;
end;

v4Candidates := [];
v4PairOrbitCounts := [];
for i in [1..Length(highClassObjects)] do
  a := Representative(highClassObjects[i]);
  C := Centralizer(K,a);
  pairOrbits := OrbitsDomain(C,eligibleDistinctPairs[i],OnPoints);
  Add(v4PairOrbitCounts,[i,Length(eligibleDistinctPairs[i]),Length(pairOrbits)]);
  for orb in pairOrbits do
    H := Group(a,orb[1]);
    if Size(H)<>4 then Error("eligible pair did not generate V4"); fi;
    Add(v4Candidates,H);
  od;
od;
v4Reps := ConjugacyDeduplicateSubgroups(v4Candidates);
v4Records := List(v4Reps,H->rec(
  order := Size(H),
  common_fixed_rank := GroupCommonFixedRank(H),
  profile := SubgroupProfile(H),
  conjugacy_class_size := Index(K,Normalizer(K,H)),
  normalizer_order := Size(Normalizer(K,H)),
  centralizer_order := Size(Centralizer(K,H)),
  generator_archived_line_images := List(GeneratorsOfGroup(H),
    x->ListPerm(TransportToArchive(x),112)),
  generator_archived_ns_matrices := List(GeneratorsOfGroup(H),
    x->ArchivedNSMatrix(TransportToArchive(x)))
));
Print("v4_pair_orbit_counts=",v4PairOrbitCounts,"\n");
Print("v4_conjugacy_class_count=",Length(v4Reps),"\n");
Print("v4_records=",List(v4Records,r->[r.common_fixed_rank,r.profile,
  r.conjugacy_class_size,r.normalizer_order,r.centralizer_order]),"\n");

ExtendHighElementary := function(inputReps)
  local candidates, H, C, hels, eligibleC, c, E;
  candidates := [];
  for H in inputReps do
    C := Centralizer(K,H);
    hels := Elements(H);
    eligibleC := Filtered(Elements(C),function(c)
      return not c in H and Order(c)=2 and
        ForAll(hels,h->HighFixedRank(h*c)<>fail);
    end);
    for c in eligibleC do
      E := Group(Concatenation(GeneratorsOfGroup(H),[c]));
      if Size(E)=2*Size(H) and GroupCommonFixedRank(E)>10 then
        Add(candidates,E);
      fi;
    od;
  od;
  return ConjugacyDeduplicateSubgroups(candidates);
end;

e8Reps := ExtendHighElementary(v4Reps);
e16Reps := ExtendHighElementary(e8Reps);
e32Reps := ExtendHighElementary(e16Reps);
Print("elementary_extension_counts=",List([e8Reps,e16Reps,e32Reps],Length),"\n");
Print("e8_records=",List(e8Reps,H->[GroupCommonFixedRank(H),SubgroupProfile(H),
  Index(K,Normalizer(K,H)),Size(Normalizer(K,H)),Size(Centralizer(K,H))]),"\n");
Print("e16_records=",List(e16Reps,H->[GroupCommonFixedRank(H),SubgroupProfile(H),
  Index(K,Normalizer(K,H)),Size(Normalizer(K,H)),Size(Centralizer(K,H))]),"\n");
Print("e32_records=",List(e32Reps,H->[GroupCommonFixedRank(H),SubgroupProfile(H),
  Index(K,Normalizer(K,H)),Size(Normalizer(K,H)),Size(Centralizer(K,H))]),"\n");

outputPath := "/private/tmp/fermat_char3_full_stabilizer_exact.g";
PrintTo(outputPath,
  "# Exact output from scripts/scratch_fermat_char3_full_stabilizer.g\n",
  "FQLineGraphIsomorphism := ",ListPerm(lineIso,112),";\n",
  "FQChamberGroupOrder := ",Size(Karchive),";\n",
  "FQProjectiveUnitarySubgroupOrder := ",Size(P),";\n",
  "FQChamberGeneratorLineImages := ",
    List(GeneratorsOfGroup(Karchive),p->ListPerm(p,112)),";\n",
  "FQChamberGeneratorNSMatrices := ",
    List(GeneratorsOfGroup(Karchive),ArchivedNSMatrix),";\n",
  "FQHighFixedRankClassRecords := ",highRecords,";\n"
  ,"FQDistinctHighPairSummary := ",pairSummary,";\n"
  ,"FQMaxDistinctPairCommonFixedRank := ",maxDistinctPairCommonFixedRank,";\n"
  ,"FQV4HighFixedRankRecords := ",v4Records,";\n"
);
Print("exact_artifact=",outputPath,"\n");

LoadPackage("json");
jsonPath := "data/enumerations/ss_fermat_quartic_char3_full_stabilizer_exact.json";
jsonPayload := rec(
  schema := "ss_fermat_quartic_char3_full_stabilizer_exact_v1",
  source_assignments := ["FQGramS","FQNSbasis","FQlineclasses","FQAF"],
  projective_unitary_subgroup_order := Size(P),
  chamber_stabilizer_order := Size(Karchive),
  conjugacy_class_count := Length(classes),
  line_graph_isomorphism := ListPerm(lineIso,112),
  chamber_generator_line_images := List(GeneratorsOfGroup(Karchive),p->ListPerm(p,112)),
  chamber_generator_ns_matrices := List(GeneratorsOfGroup(Karchive),ArchivedNSMatrix),
  high_fixed_rank_class_records := highRecords,
  v4_high_fixed_rank_records := v4Records,
  max_distinct_pair_common_fixed_rank := maxDistinctPairCommonFixedRank,
  elementary_order_8_extension_count := Length(e8Reps)
);
jsonText := GapToJsonString(jsonPayload);
jsonText := ReplacedString(jsonText,Concatenation("\\","\n"),"");
PrintTo(jsonPath,jsonText,"\n");
Print("exact_json_artifact=",jsonPath,"\n");

QUIT;
