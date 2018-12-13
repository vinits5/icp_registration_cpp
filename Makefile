TREE = KDTree/kd_tree.h KDTree/kd_tree.cpp
DATA_TYPE = KDTree/type_defs.h
ICP = ICP/ICP.h ICP/ICP.cpp
file = "icp_with_options.cpp"

all:
	make compile && make run && make results

compile:
	g++ -std=c++11 $(file) $(TREE) $(DATA_TYPE) $(ICP) -o icp_with_options

clean:
	rm icp_with_options && rm output.txt && rm transformation.txt

run:
	./icp_with_options

results:
	python display_results.py