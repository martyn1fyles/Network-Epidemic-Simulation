#Test dynamic_sbm
from Simulation.DynamicNetworks import dynamic_stochastic_block_model

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

