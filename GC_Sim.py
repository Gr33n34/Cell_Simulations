'''
Credits to 
> JTSkorik and their Monash University Winter Research Project
> Meyer-Hermann (2017), How to Simulate a Germinal Center. Methods in Molecular Biology.
'''

# Imports
import random
from enum import Enum
import pickle
import json
import csv
import os
import sys
import matplotlib.pyplot as plt

# Enumerate cell types
class CellType(Enum):
    """
    Class used to enumerate all cells in the simulation
    """
    FoBCell = 1
    Centroblast = 2
    Centrocyte = 3
    Outcell = 4

# Enumerate all cell states
class CellState(Enum):
    """
    Class used to enumerate all possible cell states
    """
    # FO B cell states
    fob_searching = 10
    fob_migrating_to_DZ = 11
    # Centroblast cell states
    cb_dividing = 12
    cb_migrating_to_LZ = 13
    # Centrocyte cell states
    cc_affinitytest = 14
    cc_back_to_DZ = 15
    cc_migrate_out = 16
    cc_apoptosis = 17

# Enumerate all cell positions
class CellLocalisation(Enum):
    '''
    Class used to enumerate all possible cell localisations
    '''
    primary_follicle = 20
    DZ = 21
    LZ = 22
    out = 23

class Params():
    """
    Class to store all relevant constant parameters of the simulation
    """
    def __init__(self):
        # Time variables
        self.dt = 0.1
        self.tmin = 0.0
        self.tmax = 200.0 #504

        # Initialisation
        self.initial_fob = 1000 #1000
        self.dist_to_DZ = 50
        self.dist_to_LZ = 50

        # Speed
        self.speed_chemokinesis = 5.0
        self.sd_chemokinesis = 0.15
        self.speed_chemotaxis = 7.5
        self.sd_chemotaxis = 0.2

        # Persistent Length time (PLT): prob of movement
        self.plt_fobcell = 0.05
        self.plt_centrocyte = 0.025
        self.plt_centroblast = 0.025

        # Probabilities (prob) and max divisions
        self.prob_divide = 0.72
        self.prob_divide_noncb = 0.0
        self.max_div_upper = 6
        self.max_div_lower = 2
        self.divide_time = 2.5
        self.begin_ag_degrad = 12
        self.end_ag_degrad = 48

        # Distance to next test
        self.dist_next_test_upper = 30
        self.dist_next_test_lower = 20
        self.dist_out = 120

        # time steps
        self.differ_to_cb_time = 0.1
        self.differ_to_cc_time = 0.1

        # Selection steps
        self.prob_fobcell_selection = 0.01 * self.dt
        self.prob_centrocyte_selection = 0.2
        self.test_delay = 0.2 # time needed for affinity testing
        self.prob_sel = self.dt * 0.05 # chance of positive signal
        self.rescue_time = 10 # time in which cell can be rescued from apoptosis

class Out():
    """
    Class to store all output of the simulation
    """
    def __init__(self,parameters):
        # Simulation time
        self.t = 0.0

        # Counter for saved data
        self.save_counter = 0

        # Available cell ids
        self.available_cell_ids = list(range(10000))

        # Lists to store cell IDs
        self.list_fob = []
        self.list_cb = []
        self.list_cc = []
        self.list_out = []

        # Lists for result plotting
        self.times = []
        self.num_fob = []
        self.num_cb = []
        self.num_cc = []
        self.num_out = []

        # Lists for result saving
        self.curr_times = []
        self.curr_num_fob = []
        self.curr_num_cb = []
        self.curr_num_cc = []
        self.curr_num_out = []

        # Cell Properties
        ## General
        self.type = {}
        self.state = {}
        self.position = {} # FO, DZ, LZ, Out
        ## FO B cell
        self.dist_to_next_test = {}
        ## Centroblasts
        self.num_divisions_to_do = {} 
        ## Centrocytes
        self.dist_to_next_test = {}
        self.clock = {} # timeout clock
        self.selfdestruc_timer = {}

# Main functions
def chemokinesis(cell_id, parameters, output):
    """
    Cell moves randomly to find next selection position
    """
    # Retrieve cell type specific movement parameters
    cell_type = output.type[cell_id]
    if cell_type == CellType.FoBCell:
        upper_speed = parameters.speed_chemokinesis * (1+parameters.sd_chemokinesis)
        lower_speed = parameters.speed_chemokinesis * (1-parameters.sd_chemokinesis) 
    elif cell_type == CellType.Centrocyte:
        upper_speed = parameters.speed_chemokinesis * (1+parameters.sd_chemokinesis)
        lower_speed = parameters.speed_chemokinesis * (1-parameters.sd_chemokinesis)
    else:
        upper_speed = None
        lower_speed = None
        print("chemokinesis: Invalid cell_type, {}".format(cell_type))

    # Move cell accordingly 
    output.dist_to_next_test[cell_id] -= parameters.dt * random.uniform(lower_speed,upper_speed)

def chemotaxis(cell_id, parameters, output):
    '''
    Moves cells between positions
    '''
    upper_speed = parameters.speed_chemotaxis * (1+parameters.sd_chemotaxis)
    lower_speed = parameters.speed_chemotaxis * (1-parameters.sd_chemotaxis)
    output.dist_to_next_test[cell_id] -= parameters.dt * random.uniform(lower_speed, upper_speed)


def testing_procedure(cell_id, parameters, output):
    """
    Function to test affinity with probability
    """
    # Test only FO B cells and centrocytes
    cell_type = output.type[cell_id]
    if cell_type == CellType.FoBCell:
        prob = parameters.prob_fobcell_selection * antigen_degradation(parameters,output)
    elif cell_type == CellType.Centrocyte:
        prob = parameters.prob_centrocyte_selection
    else:
        prob = None
        print("testing_fobcell: Invalid cell_type, {}".format(cell_type))    
    # Apply testing procedure
    if random.uniform(0,1) < prob:
        output.clock[cell_id] += parameters.test_delay
        return True
    else:
        return False

def calc_dist_for_retesting(cell_id, parameters, output):
    output.dist_to_next_test[cell_id] += random.randrange(parameters.dist_next_test_lower,parameters.dist_next_test_upper)

def antigen_degradation(parameters,output):
    if output.t < parameters.begin_ag_degrad:
        coefficient = 1
    elif output.t < parameters.end_ag_degrad:
        app = parameters.end_ag_degrad / (parameters.end_ag_degrad + parameters.begin_ag_degrad)
        coeff = (1 - app) / parameters.begin_ag_degrad
        coefficient = output.t * coeff + app
    else:
        coefficient = 0
    return coefficient

def divide(cell_id,parameters,output):
    '''
    Function to divide centroblasts
    '''
    new_cell_id = output.available_cell_ids.pop()
    output.list_cb.append(new_cell_id)
    # Creating new centroblast instance
    output.type[new_cell_id] = CellType.Centroblast
    output.state[new_cell_id] = CellState.cb_dividing
    output.position[new_cell_id] = CellLocalisation.DZ
    output.clock[new_cell_id] = output.clock[cell_id]
    output.dist_to_next_test[new_cell_id] = output.dist_to_next_test[cell_id]
    
    # Setting time for mitosis
    output.clock[new_cell_id] += parameters.divide_time
    output.clock[cell_id] += parameters.divide_time

    # Reduce number of divisions
    remaining_divisions = output.num_divisions_to_do[cell_id] - 1
    output.num_divisions_to_do[cell_id] = remaining_divisions
    output.num_divisions_to_do[new_cell_id] = remaining_divisions


def initialise_cells(parameters, output):
    '''
    Setup for simulation. Initiate cells for the simulation
    :param parameters: stores all parameters and variables in simulation
    '''
    for _ in range(parameters.initial_fob):
        cell_id = output.available_cell_ids.pop()
        output.list_fob.append(cell_id)
        output.type[cell_id] = CellType.FoBCell
        output.state[cell_id] = CellState.fob_searching
        output.position[cell_id] = CellLocalisation.primary_follicle
        output.dist_to_next_test[cell_id] = random.randrange(parameters.dist_next_test_lower,parameters.dist_next_test_upper)
        output.num_divisions_to_do[cell_id] = random.randrange(parameters.max_div_lower,parameters.max_div_upper)
        output.clock[cell_id] = 0
        output.selfdestruc_timer[cell_id] = parameters.rescue_time

def differentiate_to_cb(cell_id,parameters,output):
    '''
    Differentiate FO B cells after successful testing to centroblast
    :param cell_id: specific ID of FO B cell that will be differentiated
    :param parameters: stores all parameters of simulation
    '''
    output.type[cell_id] = CellType.Centroblast
    output.state[cell_id] = CellState.cb_dividing
    output.position[cell_id] = CellLocalisation.DZ
    output.dist_to_next_test[cell_id] = parameters.dist_to_LZ
    output.clock[cell_id] = parameters.differ_to_cb_time
    output.selfdestruc_timer[cell_id] = parameters.rescue_time

def differentiate_to_cc(cell_id,parameters,output):
    '''
    Differentiate centroblasts after end of divisions to centrocytes
    :param cell_id: specific ID of centroblast that will be differentiated
    :param parameters: stores all parameters of simulation
    '''
    output.type[cell_id] = CellType.Centrocyte
    output.state[cell_id] = CellState.cc_affinitytest
    output.position[cell_id] = CellLocalisation.LZ
    output.dist_to_next_test[cell_id] = random.randrange(parameters.dist_next_test_lower,parameters.dist_next_test_upper)
    output.clock[cell_id] = parameters.differ_to_cc_time
    output.selfdestruc_timer[cell_id] = parameters.rescue_time


def hyphasma(parameters,output,filename_output):
    '''
    Main driver function for the simulation
    :param parameters: params object that stores all parameters and variables in the simulation
    :output: current save state of the simulation
    :return:
    '''
    if output.t == 0:
        initialise_cells(parameters, output)
    
    while output.t <= parameters.tmax:
        # if 1 simulated hour has elapsed, save current state
        if output.t >= 1 * output.save_counter:
            output.save_counter += 1
            
            # saving result parameters for csv
            output.curr_times = output.t
            output.curr_num_fob = len(output.list_fob)
            output.curr_num_cb = len(output.list_cb)
            output.curr_num_cc = len(output.list_cc)
            output.curr_num_out = len(output.list_out)
            
            # saving result parameters for plotting
            output.times.append(output.t)
            output.num_fob.append(output.curr_num_fob)
            output.num_cb.append(output.curr_num_cb)
            output.num_cc.append(output.curr_num_cc)
            output.num_out.append(output.curr_num_out)

            # compiling results
            update_out_csv(output, filename_output)
            print(output.t)
            #print("Number FOB cells: {}".format(output.list_fob))
            print(str(len(output.list_fob))+' | '+str(len(output.list_cb))+' | '+str(len(output.list_cc))+' | '+str(len(output.list_out)))
        
        # iterate over FOB cells
        random.shuffle(output.list_fob)
        for i, cell_id in enumerate(output.list_fob):
            # Ignoring busy cells
            if output.clock[cell_id] > 0:
                output.clock[cell_id] -= parameters.dt
            # Testing FO B cells
            elif output.state[cell_id] == CellState.fob_searching:
                if output.dist_to_next_test[cell_id] <= 0:
                    output.clock[cell_id] += parameters.test_delay
                    if testing_procedure(cell_id, parameters, output):
                        output.state[cell_id] = CellState.fob_migrating_to_DZ
                    else:
                        calc_dist_for_retesting(cell_id, parameters, output)
                else:
                    chemokinesis(cell_id, parameters, output)
            elif output.state[cell_id] == CellState.fob_migrating_to_DZ:
                if output.dist_to_next_test[cell_id] <= 0:
                    output.list_cb.append(cell_id)
                    differentiate_to_cb(cell_id,parameters,output)
                    del (output.list_fob[i])
                else:
                    chemotaxis(cell_id, parameters, output)
            else:
                print('Error in FO B cells!')
        
        # iterate over centroblasts
        random.shuffle(output.list_cb)
        for i, cell_id in enumerate(output.list_cb):
            # Ignoring busy cells
            if output.clock[cell_id] > 0: 
                output.clock[cell_id] -= parameters.dt
            elif output.state[cell_id] == CellState.cb_dividing:
                if output.num_divisions_to_do[cell_id] > 1:
                    divide(cell_id, parameters, output)
                elif output.num_divisions_to_do[cell_id] <= 1:
                    output.state[cell_id] = CellState.cb_migrating_to_LZ
                    output.dist_to_next_test[cell_id] = parameters.dist_to_LZ
                else:
                    print('Error in dividing centroblasts!')
            elif output.state[cell_id] == CellState.cb_migrating_to_LZ:
                if output.dist_to_next_test[cell_id] > 0:
                    chemotaxis(cell_id, parameters, output)
                elif output.dist_to_next_test[cell_id] <= 0:
                    differentiate_to_cc(cell_id, parameters, output)
                    output.list_cc.append(cell_id)
                    del (output.list_cb[i])
                else:
                    print('Error in migrating centroblasts!')
            else:
                print('Error in centroblasts!')
        
        # iterate over centrocytes
        random.shuffle(output.list_cc)
        for i, cell_id in enumerate(output.list_cc):
            if output.selfdestruc_timer[cell_id] < 0:
                del (output.list_cc[i])
            elif output.state[cell_id] == CellState.cc_migrate_out:
                if output.dist_to_next_test[cell_id] < 0:
                    output.list_out.append(cell_id)
                    del (output.list_cc[i])
                else:
                    chemotaxis(cell_id, parameters, output)
            elif output.clock[cell_id] > 0:
                output.clock[cell_id] -= parameters.dt
            elif output.type[cell_id] == CellType.Centrocyte:
                if output.dist_to_next_test[cell_id] > 0:
                    chemokinesis(cell_id, parameters, output)
                elif output.dist_to_next_test[cell_id] <= 0:
                    output.clock[cell_id] += parameters.test_delay
                    if testing_procedure(cell_id,parameters,output):
                        output.state[cell_id] = CellState.cc_migrate_out  
                        output.dist_to_next_test[cell_id] = parameters.dist_out
                    else:
                        output.state[cell_id] = CellState.cc_apoptosis
                        output.selfdestruc_timer[cell_id] -= parameters.dt
                        calc_dist_for_retesting(cell_id, parameters, output)
                else:
                    print('Error in centrocyte affinity testing!')
            else:
                print('Error in centrocytes!')
            
            if output.state[cell_id] == CellState.cc_apoptosis:
                output.selfdestruc_timer[cell_id] -= parameters.dt

        # Progress time
        output.t += parameters.dt


# Helper functions
def params_to_dict(params_instance):
    '''
    Converts Params object into a dictionary with name of each variable as key for value
    :param params_instance: Params object
    :return ans: dictionary, all values of parameters
    '''
    ans = dict()
    ans['differ_to_cb_time'] = params_instance.differ_to_cb_time
    ans['differ_to_cc_time'] = params_instance.differ_to_cc_time
    ans['dist_next_test_lower'] = params_instance.dist_next_test_lower
    ans['dist_next_test_upper'] = params_instance.dist_next_test_upper
    ans['dist_to_DZ'] = params_instance.dist_to_DZ
    ans['dist_to_LZ'] = params_instance.dist_to_LZ
    ans['divide_time'] = params_instance.divide_time
    ans['dt'] = params_instance.dt
    ans['initial_fob'] = params_instance.initial_fob
    ans['max_div_lower'] = params_instance.max_div_lower
    ans['max_div_upper'] = params_instance.max_div_upper
    ans['plt_fobcell'] = params_instance.plt_fobcell
    ans['plt_centrocyte'] = params_instance.plt_centrocyte
    ans['plt_centroblast'] = params_instance.plt_centroblast
    ans['prob_centrocyte_selection'] = params_instance.prob_centrocyte_selection
    ans['prob_divide'] = params_instance.prob_divide
    ans['prob_fobcell_selection'] = params_instance.prob_fobcell_selection
    ans['prob_divide_noncb'] = params_instance.prob_divide_noncb
    ans['prob_sel'] = params_instance.prob_sel
    ans['rescue_time'] = params_instance.rescue_time
    ans['sd_chemokinesis'] = params_instance.sd_chemokinesis
    ans['sd_chemotaxis'] = params_instance.sd_chemotaxis
    ans['speed_chemokinesis'] = params_instance.speed_chemokinesis
    ans['speed_chemotaxis'] = params_instance.speed_chemotaxis
    ans['test_delay'] = params_instance.test_delay
    ans['tmin'] = params_instance.tmin
    ans['tmax'] = params_instance.tmax
    
    return ans

def dict_to_json(dictionary, filename):
    '''
    Converts dictionary object to json file and saves
    :param dictionary: dict, the dictionary intended on saving
    :param filename: str, directory to save file at
    ''' 
    with open(filename + '.json','w') as fp:
        json.dump(dictionary,fp,sort_keys=True)

def json_to_params(parameters, filename):
    '''
    Modifies a Params object to contain the exact parameters of current simulation
    :param parameters: Params object, parameters to be overwritten
    :param filename: str, directory of saved json file containing parameters
    :return:
    '''
    # open json file
    with open(filename + '.json') as fp:
        params_dict = json.load(fp)
    
    # modify current parameters object to be wanted values
    parameters['differ_to_cb_time'] = params_dict.differ_to_cb_time
    parameters['differ_to_cc_time'] = params_dict.differ_to_cc_time
    parameters['dist_next_test_lower'] = params_dict.dist_next_test_lower
    parameters['dist_next_test_upper'] = params_dict.dist_next_test_upper
    parameters['dist_to_DZ'] = params_dict.dist_to_DZ
    parameters['dist_to_LZ'] = params_dict.dist_to_LZ
    parameters['divide_time'] = params_dict.divide_time
    parameters['dt'] = params_dict.dt
    parameters['initial_fob'] = params_dict.initial_fob
    parameters['max_div_lower'] = params_dict.max_div_lower
    parameters['max_div_upper'] = params_dict.max_div_upper
    parameters['plt_fobcell'] = params_dict.plt_fobcell
    parameters['plt_centrocyte'] = params_dict.plt_centrocyte
    parameters['plt_centroblast'] = params_dict.plt_centroblast
    parameters['prob_centrocyte_selection'] = params_dict.prob_centrocyte_selection
    parameters['prob_divide'] = params_dict.prob_divide
    parameters['prob_fobcell_selection'] = params_dict.prob_fobcell_selection
    parameters['prob_divide_noncb'] = params_dict.prob_divide_noncb
    parameters['prob_sel'] = params_dict.prob_sel
    parameters['rescue_time'] = params_dict.rescue_time
    parameters['sd_chemokinesis'] = params_dict.sd_chemokinesis
    parameters['speed_chemokinesis'] = params_dict.speed_chemokinesis
    parameters['sd_chemotaxis'] = params_dict.sd_chemotaxis
    parameters['speed_chemotaxis'] = params_dict.speed_chemotaxis
    parameters['test_delay'] = params_dict.test_delay
    parameters['tmin'] = params_dict.tmin
    parameters['tmax'] = params_dict.tmax

def start_out_csv(filename):
    '''
    Creates .csv file to store output
    Writes variable name for each column
    :param filename: str, filename (and location) for .csv file
    :return:
    '''
    variable_names = ['time','num_fob','num_cb','num_cc','num_out']

    with open(filename + '.csv','w') as output_file:
        output_writer = csv.writer(output_file)
        output_writer.writerow(variable_names)

def update_out_csv(out_instance, filename):
    '''
    Writes current output values to csv file
    :param out_instance: Out object, instance of all changing variables in simulation
    :param filename: str, filename (and location) for .csv file
    :return:
    ''' 
    new_line = [out_instance.curr_times, out_instance.curr_num_fob, out_instance.curr_num_cb, out_instance.curr_num_cc, out_instance.curr_num_out]

    with open(filename + '.csv','a') as output_file:
        output_writer = csv.writer(output_file)
        output_writer.writerow(new_line)

def pickle_current_state(out_instance, simulation_name):
    """
    Saves current state of simulation using pickle for the purpose of restarting simulation.
    :param out_instance: Out object, instance of all changing variables in simulation.
    :return:
    """
    restart_data = open(simulation_name + "_restart{:04d}.pickle".format(out_instance.save_counter), "wb")
    pickle.dump(out_instance, restart_data)
    restart_data.close()


def recover_state_from_pickle(filename):
    """
    Searches through current directory to find most recent
    :return: Out object from simulation.
    """

    parameters_file = open(filename, "rb")
    output = pickle.load(parameters_file)
    return output

if __name__ == '__main__':
    '''
    Requires one input to run simulation, which is the simulation name.
    This will be used in the naming of the parameters json file, restarting pickle
    data and the output csv file.
    '''
    assert len(sys.argv) == 2, "wrong number arguments given: {}".format(len(sys.argv))

    simulation_name = sys.argv[1]

    # Find all files in current directory
    all_files = os.listdir(".")

    # Generate Params object that might be overwritten with new values
    parameters = Params()

    # Check if files exist, if not, make them:
    if simulation_name + ".json" in all_files:
        json_to_params(parameters, simulation_name)
    else:
        parameters_dict = params_to_dict(parameters)
        dict_to_json(parameters_dict, simulation_name)
    
    # Find restart files and sort so latest in final position
    restart_files = [file for file in all_files if simulation_name + "_restart" in file]
    restart_files.sort()
    if restart_files:
        # If restart files exist, load data
        output = recover_state_from_pickle(restart_files[-1])
    else:
        # otherwise, create new Out object
        output = Out(parameters)

    # If output file does not exist, create new one
    if simulation_name + ".csv" not in all_files:
        start_out_csv(simulation_name)

    # Starts/continues simulation
    hyphasma(parameters, output, simulation_name)
    
    # Plotting results
    plt.plot(output.times, output.num_cb,label='CB')
    plt.plot(output.times, output.num_cc,label='CC')
    plt.plot(output.times, output.num_out,label='Out')
    plt.legend(loc='upper right')
    plt.savefig(simulation_name+".png")