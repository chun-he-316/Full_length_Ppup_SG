# Calculate AS events
python suppa.py generateEvents -i sgd.express.gtf -o sgd.events -f ioe -e SE SS MX RI FL
python suppa.py generateEvents -i sgr.express.gtf -o sgr.events -f ioe -e SE SS MX RI FL
python suppa.py generateEvents -i ca.express.gtf -o ca.events -f ioe -e SE SS MX RI FL
