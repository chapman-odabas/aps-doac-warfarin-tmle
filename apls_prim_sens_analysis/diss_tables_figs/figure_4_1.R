library(DiagrammeR)

grViz("
digraph cohort_flow {
  graph [layout = dot, rankdir = TB]

  node [shape = box, style = rounded, fontname = Helvetica]

  aps     [label = 'APS diagnosis identified\\nN = 29,990']
  anticoag[label = 'Any anticoagulant after APS diagnosis\\nN = 11,074']
  follow  [label = 'Eligible index for full 1-year follow-up\\nN = 9,654']
  single  [label = 'Single-class initiation\\nN = 9,629']
  adult   [label = 'Adults aged >=18 years\\nN = 9,521']
  death   [label = 'Adults death_date > index date years\\nN = 9,467']
  final   [label = 'Final 1-year analytic cohort\\nN = 9,467']

  excl1 [label = 'Excluded: insufficient follow-up\\nN = 1,420']
  excl2 [label = 'Excluded: dual-class initiation\\nN = 25']
  excl3 [label = 'Excluded: age <18 years\\nN = 108']
  excl4 [label = 'Excluded: death recorded before index date\\nN = 54']

  aps -> anticoag
  anticoag -> follow
  anticoag -> excl1
  follow -> single
  follow -> excl2
  single -> adult
  single -> excl3
  adult -> excl4
  adult -> death
  death -> final
  
}
")