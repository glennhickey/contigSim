binPath=bin

all : ${binPath}/runSim.py 

clean : 
	rm -f ${binPath}/runSim.py

test: all
	python tests/allTests.py

${binPath}/runSim.py : src/runSim.py
	mkdir -p $(dir $@)
	cp src/runSim.py ${binPath}/runSim.py
	chmod +x ${binPath}/runSim.py
