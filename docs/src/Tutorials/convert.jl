using IJulia, Conda 
Conda.add("nbconvert")

folderpath = "docs/src/Tutorials"
files = readdir(folderpath)

ipynb_files = filter(file -> endswith(file, ".ipynb"), files)

for file in ipynb_files
    nbconvert = IJulia.find_jupyter_subcommand("nbconvert");
    append!(nbconvert.exec, ["--to", "markdown", "--execute", joinpath(folderpath, file) ])
    run(nbconvert)
end

