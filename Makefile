binPath=bin

all : ${binPath}/contigSim.py 

clean : 
	rm -f ${binPath}/contigSim.py

test: all
	python tests/allTests.py

${binPath}/contigSim.py : src/contigSim.py
	mkdir -p $(dir $@)
	cp src/contigSim.py ${binPath}/contigSim.py
	chmod +x ${binPath}/contigSim.py
