PROBLEMS = racetrack wetfloor

all:
	for p in $(PROBLEMS); do \
	    (cd $$p; make) \
        done

clean:
	for p in $(PROBLEMS); do \
	    (cd $$p; make clean) \
        done
	(cd engine; make)
	rm -f *~ core
