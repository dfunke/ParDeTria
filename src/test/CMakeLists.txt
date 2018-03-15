
add_executable(GridTest
	GridTest.cpp
	)

add_test(
	NAME GridTest
	COMMAND $<TARGET_FILE:GridTest>
	)
