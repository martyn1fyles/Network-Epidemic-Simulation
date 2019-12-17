#Test dynamic_sbm
from NetworkEpidemicSimulation.DynamicNetworks import dynamic_stochastic_block_model

sizes = [100, 100, 100]
probs = [[0.4, 0.001, 0.001],[0.001, 0.4, 0.001], [0.001, 0.001, 0.4]]
migration = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]
exp_par = 10
time_until = 100

test_class = dynamic_stochastic_block_model(sizes, probs, migration, exp_par, time_until)
test_class.generate_migration_times(1)

def test_generate_migration_times():

    test_class = dynamic_stochastic_block_model(sizes, probs, migration, exp_par, time_until)

    results = test_class.generate_migration_times(1)
    groups = [i[0] for i in results]
    times = [i[1] for i in results]

    #test times:
    assert times[0] == 0
    assert times[-1] > time_until

    #test groups:
    assert groups[0] == 0
    assert all([i <  4 for i in groups])

def test_get_node_current_block():
    test_class = dynamic_stochastic_block_model(sizes, probs, migration, exp_par, time_until)
    assert test_class.get_node_current_block(1) == 0
    assert test_class.get_node_current_block(101) == 1

def test_get_node_memberships():
    test_migration = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    test_class = dynamic_stochastic_block_model(sizes, probs, test_migration, exp_par, time_until)
    assert all([membership == 0 for membership in test_class.get_node_memberships(1)])

def test_get_node_next_migration_time():
    # Test was written before we supported non-exponential waiting times
    # Only test that the vector is of approximately the right form.
    test_class = dynamic_stochastic_block_model(sizes, probs, migration, exp_par, time_until)
    assert len(test_class.get_next_migration_times()) == 300
    assert all([ 0 < time < 100 for time in test_class.get_next_migration_times()])

def test_determine_nodes_to_migrate():
    test_class = dynamic_stochastic_block_model(sizes, probs, migration, exp_par, time_until)
    
    # We overwrite the next migration time variable for every node in the class with 5
    [test_class.G.nodes[node].update({"Next Migration Time": 5}) for node in test_class.G.nodes()]

    assert test_class.determine_nodes_to_migrate(4) == []
    assert test_class.determine_nodes_to_migrate(6) == list(test_class.G.nodes())
    

def test_perform_migration_event():
    # We set up the migration probabilities so that we know prob 1 which block the node will migrate to
    test_migration = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    test_class = dynamic_stochastic_block_model(sizes, probs, test_migration, exp_par, time_until)
    
    # We overwrite the next migration time variable for every node in the class with 5
    [test_class.G.nodes[node].update({"Next Migration Time": 5}) for node in test_class.G.nodes()]

    test_class.perform_migration_event(1)
    test_class.perform_migration_event(101)
    test_class.perform_migration_event(201)

    # Nodes in block 0 move to block 1, block 1 -> block 2, block 2 -> block 0
    assert test_class.get_node_current_block(1) == 1
    assert test_class.get_node_current_block(101) == 2
    assert test_class.get_node_current_block(201) == 0

def test_update_edges():
    # update_edges is called during perform_migration_event

    # Network has 3 fully connected subgraphs that are disconnected
    test_probs = [[1, 0, 0],[0, 1, 0], [0, 0, 1]]
    test_migration = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    test_class = dynamic_stochastic_block_model(sizes, test_probs, test_migration, exp_par, time_until)

    # We overwrite the next migration time variable for every node in the class with 5
    [test_class.G.nodes[node].update({"Next Migration Time": 5}) for node in test_class.G.nodes()]

    # We perform the migration event for node 201, it will be moving to block 0, and should be connected to nodes 0:100
    test_class.perform_migration_event(201)

    assert list(test_class.G.neighbors(201)) == list(range(100))

def test_increment_network():
        # Network has 3 fully connected subgraphs that are disconnected
    test_probs = [[1, 0, 0],[0, 1, 0], [0, 0, 1]]
    test_migration = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    test_class = dynamic_stochastic_block_model(sizes, test_probs, test_migration, 50, time_until)

    # We overwrite the next migration time variable for every node in the class with 5
    [test_class.G.nodes[node].update({"Next Migration Time": 5}) for node in test_class.G.nodes()]
    # We overwrite the next migration time of 201 to be 0.1
    test_class.G.nodes[201].update({"Next Migration Time": 0.1})
    # The network is incremented forwards by a time of 
    test_class.increment_network(0.2)

    # We perform the migration event for node 201, it will be moving to block 0, and should be connected to nodes 0:100
    # It possible that the node could be migrated twice but really unlikely
    assert list(test_class.G.neighbors(201)) == list(range(100))
