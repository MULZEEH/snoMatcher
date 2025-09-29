awk '{print "rule " $0 " :\n\t input: \n\t output: \n\t script: \n\n#\n#\n#\n\n "}' todo.txt >> Snakefile && cat todo.txt | xargs -I {} touch scripts/"{}.R"
