library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

flowchart_lr <- DiagrammeR::grViz("digraph {

# initiate graph
graph [layout = dot, rankdir = LR, label = '',labelloc = t]

# global node settings
node [shape = rectangle, style = filled, fillcolor = Bisque, fontname = 'Helvetica']

# label nodes
data1 [label = 'GISAID\n US genome sequences \n metadata', shape = folder, fillcolor = PaleGreen]
data2 [label = 'covidestim\n State infections \n estimates', shape = folder, fillcolor = PaleGreen]
process [label =  'Omicron* \n Variants. infections \n estimates']
statistical [label = 'effective reproduction number (Rt) \n estimates']
results [label= 'Results \n - Attack Rate \n - Incidence \n - (Rt)']
metaanalysis [label = 'Meta-analysis \n Explanatory variables to variability on \n (Rt), A.R., Inc., etc.']

# edge definitions with the node IDs
{data1 data2}  -> process -> statistical
statistical -> {results metaanalysis}
metaanalysis -> results
}")

flowchart_lr |> 
  export_svg() |> 
  charToRaw() |> 
  rsvg_png("Output/Plots/flowchart_lr.png")

flowchart_lr |> 
  export_svg() |> 
  charToRaw() |> 
  rsvg_pdf("~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/flowchart_lr.pdf")

flowchart_tb <- DiagrammeR::grViz("digraph {

# initiate graph
graph [layout = dot, rankdir = TB, label = '',labelloc = t]

# global node settings
node [shape = rectangle, style = filled, fillcolor = Bisque, fontname = 'Helvetica']

# label nodes
data1 [label = 'GISAID\n US genome sequences \n metadata', shape = folder, fillcolor = PaleGreen]
data2 [label = 'covidestim\n State infections \n estimates', shape = folder, fillcolor = PaleGreen]
process [label =  'Omicron* \n Variants infections \n estimates']
statistical [label = 'effective reproduction number (Rt) \n estimates']
results [label= 'Results \n - Attack Rate \n - Incidence \n - (Rt)']
metaanalysis [label = 'Meta-analysis \n on (Rt) values']

# edge definitions with the node IDs
{data1 data2}  -> process -> statistical
statistical -> {results metaanalysis}
metaanalysis -> results
}")

flowchart_tb |> 
  export_svg() |> 
  charToRaw() |> 
  rsvg_png("Output/Plots/flowchart_tb.png")

flowchart_tb |> 
  export_svg() |> 
  charToRaw() |> 
  rsvg_pdf("~/Dropbox/GLab_team/papers/2023_Omicron-infections/Figures/flowchart_tb.pdf")
