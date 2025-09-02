
params = {'text.usetex': True,
          'text.latex.preamble': r'\usepackage{amsfonts}',
          'font.family': 'serif',
          'font.sans-serif':'Computer Modern Roman'}


fontsizeLabel=20
fontsizeLeg=18
fontsizeAnnotate=12
fontsizeSlope=16


names = { "dSh_Quad_DG0_Godunov"       : {"mesh":"Quads","EF":r"$\mathbf{RT}^1$","flux":"Explicit","marker":'D',"markerfacecolor":None,"markersize":9,"ls":'dotted',"color":'blue'},
 
          "dSh_Quad_DG0_LaxFriedrich"  : {"mesh":"Quads","EF":r"$\mathbf{RT}^1$","flux":"Implicit","marker":'v',"markerfacecolor":None,"markersize":9,"ls":'solid',"color":'red'},
          "dSh_Quad_DG1_LaxFriedrich"  : {"mesh":"Quads","EF":r"$\mathbb{P}_1$","flux":"Rusanov","marker":'>',"markerfacecolor":'None',"markersize":9,"ls":'dashed',"color":'orange'},
          "dSh_Quad_DG2_LaxFriedrich"  : {"mesh":"Quads","EF":r"$\mathbb{P}_2$","flux":"Rusanov","marker":'D',"markerfacecolor":'None',"markersize":9,"ls":'dashed',"color":'magenta'},
          "dQ_Quad_DG0_Godunov"        : {"mesh":"Quads","EF":r"$\left(\mathrm{d}\mathbb{Q}_0\right)^2$","flux":"Godunov","marker":'s',"markerfacecolor":'None',"markersize":9,"ls":'dotted',"color":'black'},
          "dQ_Quad_DG0_LaxFriedrich"   : {"mesh":"Quads","EF":r"$\mathbf{RT}^1$","flux":"Explicit Hodge-Laplacian","marker":'v',"markerfacecolor":None,"markersize":9,"ls":'solid',"color":'red'}}
