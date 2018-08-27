cmake_minimum_required(VERSION 2.6)
project(pixfit)

include_directories(${VIF_INCLUDE_DIRS})


add_executable(pixfit-extract pixfit-extract.cpp)
target_link_libraries(pixfit-extract ${VIF_LIBRARIES})
install(TARGETS pixfit-extract DESTINATION bin)

add_executable(pixfit-cutout pixfit-cutout.cpp)
target_link_libraries(pixfit-cutout ${VIF_LIBRARIES})
install(TARGETS pixfit-cutout DESTINATION bin)

add_executable(pixfit-findprior pixfit-findprior.cpp)
target_link_libraries(pixfit-findprior ${VIF_LIBRARIES})
install(TARGETS pixfit-findprior DESTINATION bin)

add_executable(pixfit-addflx pixfit-addflx.cpp)
target_link_libraries(pixfit-addflx ${VIF_LIBRARIES})
install(TARGETS pixfit-addflx DESTINATION bin)

add_executable(pixfit-addtdcst pixfit-addtdcst.cpp)
target_link_libraries(pixfit-addtdcst ${VIF_LIBRARIES})
install(TARGETS pixfit-addtdcst DESTINATION bin)

add_executable(pixfit-combine pixfit-combine.cpp)
target_link_libraries(pixfit-combine ${VIF_LIBRARIES})
install(TARGETS pixfit-combine DESTINATION bin)

add_executable(pixfit-gfit pixfit-gfit.cpp)
target_link_libraries(pixfit-gfit ${VIF_LIBRARIES})
install(TARGETS pixfit-gfit DESTINATION bin)

add_executable(pixfit-prep pixfit-prep.cpp)
target_link_libraries(pixfit-prep ${VIF_LIBRARIES})
install(TARGETS pixfit-prep DESTINATION bin)

add_executable(pixfit-reg2fits pixfit-reg2fits.cpp)
target_link_libraries(pixfit-reg2fits ${VIF_LIBRARIES})
install(TARGETS pixfit-reg2fits DESTINATION bin)

add_executable(pixfit-residual pixfit-residual.cpp)
target_link_libraries(pixfit-residual ${VIF_LIBRARIES})
install(TARGETS pixfit-residual DESTINATION bin)

add_executable(pixfit-show pixfit-show.cpp)
target_link_libraries(pixfit-show ${VIF_LIBRARIES})
install(TARGETS pixfit-show DESTINATION bin)
