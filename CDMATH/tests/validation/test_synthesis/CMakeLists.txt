
SET(NOTEBOOK
  ./test_synthesis_nbook.ipynb
  )
SET(TABLES
  ./test_synthesis_all.csv
  ./test_synthesis_short.csv
  ./test_synthesis_all.xlsx
  ./test_synthesis_short.xlsx
  )

file(COPY ${NOTEBOOK} ${TABLES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


