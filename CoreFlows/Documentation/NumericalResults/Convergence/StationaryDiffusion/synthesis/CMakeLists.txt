
SET(NOTEBOOK
  ./convergence_synthesis_StationaryDiffusion.ipynb
  )
SET(TABLES
  ./convergence_synthesis_StationaryDiffusion_all.csv
  ./convergence_synthesis_StationaryDiffusion_short.csv
  ./convergence_synthesis_StationaryDiffusion_all.xlsx
  ./convergence_synthesis_StationaryDiffusion_short.xlsx
  )

file(COPY ${NOTEBOOK} ${TABLES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target ( notebook COMMAND jupyter-nbconvert --to notebook --execute convergence_synthesis_StationaryDiffusion.ipynb )
