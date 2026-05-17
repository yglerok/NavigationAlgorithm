cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set(PROJECT_NAME NavigationAlgorithm)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

find_package(Gnuplot REQUIRED)

# Создание исполняемого файла
add_executable(${PROJECT_NAME} src/main.cpp)

if(TARGET Gnuplot::gnuplot)
    target_link_libraries(${PROJECT_NAME} PRIVATE Gnuplot::gnuplot)
else()
    # Резервный вариант: подключаем заголовочные файлы и библиотеку напрямую
    target_include_directories(${PROJECT_NAME} PRIVATE ${Gnuplot_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${Gnuplot_LIBRARIES})
endif()

# Передаём путь к исполняемому файлу Gnuplot как макрос для использования в коде
target_compile_definitions(${PROJECT_NAME} PRIVATE
    GNUPLOT_EXECUTABLE="${Gnuplot_EXECUTABLE}"
)