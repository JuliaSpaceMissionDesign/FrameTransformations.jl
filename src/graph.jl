# common export
export find_path, connect!, register!

"""
    find_path(g::G, from, to)
    find_path(from, to)  # module specific

Find shortest path in a graph `g` that starts at `from` and ends at `to`.

!!!! warning 
    This method is abstract!
"""
function find_path end

"""
    connect!(g::G, n1, n2)
    connect!(n1, n2)  # module specific

Connect a node `n2` to the node `n1` in the graph `g`.

!!!! warning 
    This method is abstract!
"""
function connect! end

"""
    register!(g::G, n)
    register!(n)  # module specific 
    
Register a new node in the graph `g`.

!!!! warning 
    This method is abstract!
"""
function register! end
