
# copy it as it

file(GLOB MY_FILES
  "./*.png"
)
file(COPY ${MY_FILES} DESTINATION .) # BUILD
file(COPY ${MY_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/images) # INSTALL