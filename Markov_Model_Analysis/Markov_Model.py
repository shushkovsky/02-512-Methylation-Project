import pandas as pd 
import numpy as np
import copy
import matplotlib.pyplot as plt
import scipy
from scipy import stats 
from scipy.stats import linregress 

def construct_state_matrix(df,methyl_names,age,methyl,locations):
    states_matrix_temp = []

    #determine thresholds for age buckets based on quartile analysis
    age_buckets = []
    for buckets in range(age):
        age_buckets += [(df.age.quantile(1/age*buckets),df.age.quantile(1/age*(buckets+1)))]

    #come up with datastructure for state matrix
    for elems in range(len(age_buckets)):
        states_matrix_mini = []
        for methylgroups in range(methyl):
            states_matrix_mini += [[]]
        states_matrix_temp += [states_matrix_mini]
        states_matrix = copy.deepcopy(states_matrix_temp)

    #add rows of the dataframe based on what age bucket they are in
    for index, row in df.iterrows():
        states_matrix_index = 0
        for (lb,ub) in age_buckets:
            if row['age'] <= ub and row['age'] >= lb:
                for methylsites in range(methyl):
                    states_matrix_temp[states_matrix_index][methylsites] += [row[methyl_names[methylsites]]]
            states_matrix_index += 1

    #modify the methylation values so they are partitioned into 'number of locations (10)' equally distanced points that reflect the entire range of datasets' methylation sites of that age group
    for agegroup_index in range(len(states_matrix_temp)):
        for methylgroup_index in range(len(states_matrix_temp[agegroup_index])):
            min_methylvalue = min(states_matrix_temp[agegroup_index][methylgroup_index])
            max_methylvalue = max(states_matrix_temp[agegroup_index][methylgroup_index])
            step = (max_methylvalue-min_methylvalue)/locations
            modified_list_methylvalues = []
            for positions in range(locations):
                modified_list_methylvalues += [min_methylvalue + step*positions]
            states_matrix[agegroup_index][methylgroup_index] += modified_list_methylvalues

    return states_matrix, age_buckets

def construct_state_matrix_prune(df,methyl_names,age,methyl,locations):
    #determine thresholds for age buckets based on quartile analysis
    age_buckets = []
    for buckets in range(age):
        age_buckets += [(df.age.quantile(1/age*buckets),df.age.quantile(1/age*(buckets+1)))]
    
    #compute the buckets for each methylation site (note it doesn't change as we change our current age buckets)
    methyl_buckets = []
    for methyl_buckets_index in range(methyl):
        methyl_name_at_index = methyl_names[methyl_buckets_index]
        methyl_buckets_mini = []
        for location_index in range(locations):
            methyl_buckets_mini += [(getattr(df,methyl_name_at_index).quantile(1/locations*location_index),getattr(df,methyl_name_at_index).quantile(1/locations*(location_index+1)))]
        methyl_buckets += [methyl_buckets_mini]
    
    #come up with datastructure for state matrix_temp
    states_matrix_temp = []
    for elems in range(len(age_buckets)):
        states_matrix_small = []
        for methylgroup in range(methyl):
            states_matrix_mini = []
            for methylgroup_locations in range(locations):
                states_matrix_mini += [[]]
            states_matrix_small += [states_matrix_mini]
        states_matrix_temp += [states_matrix_small]
    states_matrix = []
    for age_elems in range(len(age_buckets)):
        states_matrix_mini = []
        for methyl_groups in range(methyl):
            states_matrix_mini += [[]]
        states_matrix += [states_matrix_mini]

    print(len(states_matrix),len(states_matrix[0]),len(states_matrix[0][0]))
    #add rows of the dataframe based on what age bucket they are in
    for index, row in df.iterrows():
        states_matrix_age_index = 0
        for (lb,ub) in age_buckets:
            if row['age'] <= ub and row['age'] >= lb:
                for methylsites in range(methyl):
                    states_matrix_methyl_index = 0
                    for (mlb,mub) in methyl_buckets[methylsites]:
                        if row[methyl_names[methylsites]] <= mub and row[methyl_names[methylsites]] >= mlb: 
                            states_matrix_temp[states_matrix_age_index][methylsites][states_matrix_methyl_index] += [row[methyl_names[methylsites]]]
                        states_matrix_methyl_index += 1
            states_matrix_age_index += 1
    
    #For each age group compute the average value for methylation site which determines what our states will have (dim: age x methyl x locations)
    for age_group in range(age):
        for methyl_site in range(methyl):
            averaged_values = []
            for methyl_positions in range(locations):
                if len(states_matrix_temp[age_group][methyl_site][methyl_positions]) != 0:
                    average_val = sum(states_matrix_temp[age_group][methyl_site][methyl_positions])/len(states_matrix_temp[age_group][methyl_site][methyl_positions])
                else: 
                    average_val = 0
                averaged_values += [average_val]
            states_matrix[age_group][methyl_site] = averaged_values

    return states_matrix, age_buckets, methyl_buckets

def construct_state_matrix_prune_2(df,methyl_names,age,methyl,locations):
    #determine thresholds for age buckets - equally distributed
    age_buckets = [(0,12),(12,24),(24,36),(36,48),(48,60),(60,72),(72,84),(84,96),(96,108),(108,120)]
    
    #compute the buckets for each methylation site - new method
    methyl_buckets = []
    for methyl_buckets_index in range(methyl):
        methyl_name_at_index = methyl_names[methyl_buckets_index]
        methyl_buckets_mini = []
        for location_index in range(locations):
            methyl_buckets_mini += [(8/locations * location_index - 4, 8/locations * (location_index+1) - 4)]
        methyl_buckets += [methyl_buckets_mini]

    #come up with datastructure for state matrix_temp
    states_matrix_temp = []
    for elems in range(len(age_buckets)):
        states_matrix_small = []
        for methylgroup in range(methyl):
            states_matrix_mini = []
            for methylgroup_locations in range(locations):
                states_matrix_mini += [[]]
            states_matrix_small += [states_matrix_mini]
        states_matrix_temp += [states_matrix_small]
    states_matrix = []
    for age_elems in range(len(age_buckets)):
        states_matrix_mini = []
        for methyl_groups in range(methyl):
            states_matrix_mini += [[]]
        states_matrix += [states_matrix_mini]
    print(len(states_matrix),len(states_matrix[0]),len(states_matrix[0][0]))

    #add rows of the dataframe based on what age bucket they are in
    for index, row in df.iterrows():
        states_matrix_age_index = 0
        for (lb,ub) in age_buckets:
            if row['age'] <= ub and row['age'] >= lb:
                for methylsites in range(methyl):
                    states_matrix_methyl_index = 0
                    for (mlb,mub) in methyl_buckets[methylsites]:
                        if row[methyl_names[methylsites]] <= mub and row[methyl_names[methylsites]] >= mlb: 
                            states_matrix_temp[states_matrix_age_index][methylsites][states_matrix_methyl_index] += [row[methyl_names[methylsites]]]
                        states_matrix_methyl_index += 1
            states_matrix_age_index += 1
    
    #For each age group compute the average value for methylation site which determines what our states will have (dim: age x methyl x locations)
    for age_group in range(age):
        for methyl_site in range(methyl):
            averaged_values = []
            for methyl_positions in range(locations):
                if len(states_matrix_temp[age_group][methyl_site][methyl_positions]) != 0:
                    average_val = sum(states_matrix_temp[age_group][methyl_site][methyl_positions])/len(states_matrix_temp[age_group][methyl_site][methyl_positions])
                else: 
                    average_val = 0
                averaged_values += [average_val]
            states_matrix[age_group][methyl_site] = averaged_values

    return states_matrix, age_buckets, methyl_buckets


def construct_transition_matrix(state_matrix, age, methyl, locations):
    transition_matrix = []

    #come up with datastructure for transition matrix
    for ages in range(age-1):
        age_list = []
        for methylsites in range(locations):
            methylist = []
            for positions in range(locations):
                methylist += []
            age_list += [methylist]
        transition_matrix += [age_list]

    #find the incoming probabilities of transition for each of the nodes and normalize it so that all outgoing edges from one node's probability is 1 (dim: age-1 x locationA x methyl x locationB)
    index = 0 
    for ages_2 in range(age-1):
        for methylsites_2 in range(methyl):
            for positions_2 in range(locations):
                value_at_index = []
                for positions_2_newage in range(locations):
                    sum_of_transition = 0
                    for positions_2_newage_2 in range(locations):
                        if state_matrix[ages_2+1][methylsites_2][positions_2_newage_2] != 0:
                            sum_of_transition += (state_matrix[ages_2][methylsites_2][positions_2_newage]/state_matrix[ages_2+1][methylsites_2][positions_2_newage_2])
                        else: 
                            sum_of_transition += 0
                        index += 1
                    
                    if sum_of_transition != 0 and state_matrix[ages_2+1][methylsites_2][positions_2_newage] != 0:
                        value_at_index += [1/sum_of_transition*state_matrix[ages_2][methylsites_2][positions_2_newage]/state_matrix[ages_2+1][methylsites_2][positions_2_newage]]
                    else:
                        value_at_index += [0]
                transition_matrix[ages_2][positions_2] += [value_at_index]
    
    #modify it so we can just average the methylation values to one transition value (dim age-1 x locationA x locationB)
    transition_matrix_mod = []
    for age_group in transition_matrix:
        transition_matrix_mod_sublist = []
        for outgoing_locations in age_group:
            transition_matrix_mod_outgoing_list = []
            transposed_list = np.array(outgoing_locations).T.tolist()
            for index in range(len(transposed_list)):
                '''
                sum_of_squares = 0
                for elem in transposed_list[index]:
                    sum_of_squares += (elem**2)
                transition_matrix_mod_outgoing_list += [sum_of_squares/len(transposed_list[index])]
                '''
                transition_matrix_mod_outgoing_list += [sum(transposed_list[index])/len(transposed_list[index])]
            transition_matrix_mod_sublist += [transition_matrix_mod_outgoing_list]
        transition_matrix_mod += [transition_matrix_mod_sublist]
    return transition_matrix_mod


def list_of_predicted_actual_age(df,state_matrix,age_groups,methylsites,age_buckets):
    #transpose the state_matrix
    state_matrix_transposed = []
    for agegroup in state_matrix:
        state_matrix_transposed += [np.array(agegroup).T.tolist()]

    #get list of predicted ages and actual age column of df
    list_of_predicted_ages = []
    list_of_actual_ages = []
    for index, row in df.iterrows():
        list_of_actual_ages += [row['age']]
        predicted_age = 0
        min_distance = float('inf')
        for age_index in range(len(state_matrix_transposed)):
            for location_index in range(len(state_matrix_transposed[age_index])):
                distance = 0
                for methyl_index in range(len(methylsites)):
                    #add distances by absolute value distance from actual methylation site to approximated methylation site value
                    distance += abs(row[methylsites[methyl_index]] - state_matrix_transposed[age_index][location_index][methyl_index])
                #smallest distance away from actual value
                if distance < min_distance:
                    min_distance = distance
                    age_min,age_max = age_buckets[age_index]
                    predicted_age = (age_max+age_min)/2 
        list_of_predicted_ages += [predicted_age]
    
    return list_of_predicted_ages,list_of_actual_ages


def graph_ages(predicted_ages,actual_ages,condition):
    #total number of people
    enumerated_people = []
    for i in range(len(predicted_ages)):
        enumerated_people += [i+1]

    title = ""
    if condition == False:
        title = "Healthy Patients Prediction and Actual Ages"
    if condition == True:
        title = "Diabetic Patients Prediction and Actual Ages"

    plt.title(title)
    plt.plot(enumerated_people,predicted_ages, label = 'Predicted Ages')
    plt.plot(enumerated_people,actual_ages, label = 'Actual Ages')
    plt.xlabel('Patient Number')
    plt.legend()
    plt.show()

    plt.clf()
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(predicted_ages, actual_ages)
    print("Age Prediction Best Fit Line Parameters")
    print('Slope:', slope, 'Intercept:', intercept, 'R-Value:', r_value, 'P-Value:', p_value,'Standard Error:', std_err)
    plt.plot(actual_ages, predicted_ages, 'x', label = 'Prediction over Actual Ages')
    best_fit_line_points = [slope*i+intercept for i in actual_ages]
    plt.plot(actual_ages, best_fit_line_points, label = 'Linear Regression')
    best_line = [1*i+0 for i in actual_ages]
    plt.plot(actual_ages, best_line, color = 'r', label = 'y = x')
    plt.title(title)
    plt.errorbar(actual_ages, predicted_ages, yerr= std_err, ls = 'none')
    plt.xlabel('Actual Age')
    plt.ylabel('Predicted Age')
    plt.legend()
    plt.show()
    return

def average_distance_apart(predicted_age_list, actual_age_list):
    distance_array = []
    for index in range(len(predicted_age_list)):
        distance_array += [abs(actual_age_list[index]-predicted_age_list[index])]
    return sum(distance_array)/len(distance_array)

def markov_chain(state_matrix, transition_matrix):
    markov_chain_states = []
    markov_chain_transitions_previous = []
    markov_chain_transitions = []

    for age_index in range(len(state_matrix)):
        state_value = 0
        for methylation_index in range(len(state_matrix[age_index])):
            for location_index in range(len(state_matrix[age_index][methylation_index])):
                state_value += state_matrix[age_index][methylation_index][location_index]
        markov_chain_states += [state_value/(len(state_matrix[0][0])*len(state_matrix[0]))]
    
    for age_index_2 in range(len(transition_matrix)):
        transition_value = 0
        for location_1 in range(len(transition_matrix[age_index_2])):
            for location_2 in range(len(transition_matrix[age_index_2][location_1])):
                transition_value += abs(transition_matrix[age_index_2][location_1][location_2])
        markov_chain_transitions_previous += [transition_value/len(state_matrix[0][0])]

    for elem in markov_chain_transitions_previous:
        markov_chain_transitions += [elem/sum(markov_chain_transitions_previous)]

    return [markov_chain_states, markov_chain_transitions]

def condition_prediction(state_matrix_healthy, transition_matrix_healthy, state_matrix_disease, transition_matrix_disease, age_buckets_healthy, methyl_buckets_healthy, age_buckets_diseased, methyl_buckets_diseased, past_methylation_data, past_age, past_condition, current_methylation_data, current_age):
    #determine where to start from with which matrices
    patient_health = past_condition
    if patient_health == False:
        state_matrix = state_matrix_healthy
        transition_matrix = transition_matrix_healthy
        age_buckets = age_buckets_healthy
        methyl_buckets = methyl_buckets_healthy
    if patient_health == True:
        state_matrix = state_matrix_disease
        transition_matrix = transition_matrix_disease
        age_buckets = age_buckets_diseased
        methyl_buckets = methyl_buckets_diseased

    #compute error terms whether healthy or have disease
    error_of_healthy = error_of_prediction_actual(past_age, current_age, methyl_buckets_healthy, age_buckets_healthy, state_matrix_healthy, transition_matrix_healthy, past_methylation_data, current_methylation_data)
    error_of_disease = error_of_prediction_actual(past_age, current_age, methyl_buckets_diseased, age_buckets_diseased, state_matrix_disease, transition_matrix_disease, past_methylation_data, current_methylation_data)

    #determine condition of individual
    patient_condition = 0
    if abs(error_of_healthy) > abs(error_of_disease): 
        patient_condition = -1
    if abs(error_of_disease) > abs(error_of_healthy):
        patient_condition = 1
    return patient_condition

def error_of_prediction_actual(past_age, current_age, methyl_buckets, age_buckets, state_matrix, transition_matrix, past_methylation_data, current_methylation_data):
    #determine our age_index to start from
    age_index_past = -1
    for age_index in range(len(age_buckets)):
        age_lb, age_ub = age_buckets[age_index]
        if past_age > age_lb and past_age < age_ub:
            age_index_past = age_index

    #determine our methylation index to start from 
    sum_at_methylation_site = [0]*len(methyl_buckets[0])
    for methyl_index in range(len(methyl_buckets)):
        for location_index in range(len(methyl_buckets[methyl_index])):
            methyl_lb, methyl_ub = methyl_buckets[methyl_index][location_index]
            if past_methylation_data[methyl_index] > methyl_lb and past_methylation_data[methyl_index] < methyl_ub:
                sum_at_methylation_site[location_index] += 1
    location_index_past = -1
    location_score = -1
    for sum_methylation_index in range(len(sum_at_methylation_site)):
        if sum_at_methylation_site[sum_methylation_index] > location_score:
            location_score = sum_at_methylation_site[sum_methylation_index]
            location_index_past = sum_methylation_index

    #indices_traveled = [location_index_past]
    loop_location_index = location_index_past
    loop_age_index = age_index_past
    if age_buckets[-1][1] < current_age or age_buckets[0][0] > current_age:
        return float('inf')
    while current_age < age_buckets[loop_age_index][0] or current_age > age_buckets[loop_age_index][1]:
        next_location_index = transition_matrix[loop_age_index][loop_location_index].index(max(transition_matrix[loop_age_index][loop_location_index]))
        #indices_traveled += [next_location_index]
        loop_location_index = next_location_index 
        loop_age_index += 1

    #error score computed 
    score = 0
    for methyl_site in range(len(current_methylation_data)):
        score += (current_methylation_data[methyl_site] - state_matrix[loop_age_index][methyl_site][loop_location_index])

    return score


def predictions_of_conditions(longitudinal_testing_df, state_matrix_healthy, transition_matrix_healthy, state_matrix_diseased, transition_matrix_diseased, age_buckets_healthy, methyl_buckets_healthy, age_buckets_diseased, methyl_buckets_diseased):
    #convert dataframe to list
    numpy_array = longitudinal_testing_df.to_numpy()
    list_of_patients = []
    for elem in numpy_array:
        list_of_patients += [list(elem)]

    #list of conditions
    predicted_condition = []
    actual_condition = []   

    #iterate by only locating the previously seen version of that patient
    seen_patients = []
    past_methylation_data = []
    past_age = -1
    past_condition = False 
    current_methylation_data = []
    current_age = -1
    current_condition = False 
    for elems in list_of_patients:
        if elems not in seen_patients:
            seen_patients += [elems]
            past_methylation_data = elems[2:-1]
            past_age = elems[1]
            past_condition = elems[-1]
        if elems in seen_patients:
            current_methylation_data = elems[2:-1]
            current_age = elems[1]
            current_condition = elems[-1]
            #actual conditions of this patient
            if current_condition == True:
                actual_condition += [-1]
            if current_condition == False:
                actual_condition += [1]
            #predicted conditions of this patient
            predicted_condition += [condition_prediction(state_matrix_healthy, transition_matrix_healthy, state_matrix_diseased, transition_matrix_diseased, age_buckets_healthy, methyl_buckets_healthy, age_buckets_diseased, methyl_buckets_diseased, past_methylation_data, past_age, past_condition, current_methylation_data, current_age)]
            past_methylation_data = []
            past_age = -1
            past_condition = False
            current_methylation_data = []
            current_age = -1
            current_condition = False
    return predicted_condition, actual_condition

def graph_conditions(predicted_condition, actual_condition):
    x_coordinates = []
    for elems in range(len(predicted_condition)):
        x_coordinates += [elems]
    
    plt.plot(x_coordinates, predicted_condition, label = "Predicted Conditions")
    plt.plot(x_coordinates, actual_condition, label = "Actual Conditions") 
    plt.legend()
    plt.show()

    plt.clf()
    barWidth = 0.15
    bar_plot_predicted = [predicted_condition.count(1), predicted_condition.count(-1)]
    bar_plot_actual = [actual_condition.count(1), actual_condition.count(-1)]
    bar_plot_x_1 = [0,1]
    bar_plot_x_2 = [x+barWidth for x in bar_plot_x_1]
    plt.bar(bar_plot_x_1, bar_plot_predicted, width = barWidth, label = 'Predicted')
    plt.bar(bar_plot_x_2, bar_plot_actual, width = barWidth, label = 'Actual')
    plt.xticks([r+barWidth for r in range(len(bar_plot_actual))],['Healthy Patients', 'Diabetic Patients'])
    plt.legend()
    plt.show()
    print('Healthy:', 'Predicted', predicted_condition.count(1), 'Actual', actual_condition.count(1))
    print('Diabetes:', 'Predicted', predicted_condition.count(-1), 'Actual', actual_condition.count(-1))
    return 

def main():
    #modify to wherever you have kept the csv file 
    df_disease = pd.read_csv("/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Input_Datasets/Diabetes_Set.csv")
    df_healthy = pd.read_csv("/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Input_Datasets/Healthy_Set.csv")
    list_of_methylation_sites = list(df_healthy.columns)
    list_of_methylation_sites = list_of_methylation_sites[1:len(list_of_methylation_sites)]
    number_of_methylation_sites = len(list_of_methylation_sites)
    number_of_age_groups = 10
    number_of_locations_at_each_age_group = 100

    print('Dimensions:')
    print('Number of Methylation Sites', number_of_methylation_sites)
    print('Number of Age Groups', number_of_age_groups)
    print('Number of Location Positions at each Age Group', number_of_locations_at_each_age_group)

    #Predict Healthy Age Pre-Prune
    print('\n Pre-Prune-Healthy')
    #make markov model
    state_matrix_healthy_pre_prune, age_buckets_healthy_pre_prune = construct_state_matrix(df_healthy,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_healthy_pre_prune = construct_transition_matrix(state_matrix_healthy_pre_prune, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_healthy_pre_prune = pd.DataFrame(state_matrix_healthy_pre_prune)
    df_state_healthy_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_healthy_pre_prune.csv')
    df_transition_healthy_pre_prune = pd.DataFrame(transition_matrix_healthy_pre_prune)
    df_transition_healthy_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_healthy_pre_prune.csv')
    #predict age 
    predicted_ages_healthy_pre_prune, actual_ages_healthy_pre_prune = list_of_predicted_actual_age(df_healthy,state_matrix_healthy_pre_prune,number_of_age_groups,list_of_methylation_sites,age_buckets_healthy_pre_prune)
    print('Healthy_Pre_Prune_Error', average_distance_apart(predicted_ages_healthy_pre_prune, actual_ages_healthy_pre_prune))
    graph_ages(predicted_ages_healthy_pre_prune, actual_ages_healthy_pre_prune,False)

    #Predict Diseased Age Pre-Prune
    print('\n Pre-Prune-Disease')
    #make markov model
    state_matrix_diseased_pre_prune, age_buckets_diseased_pre_prune = construct_state_matrix(df_disease,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_diseased_pre_prune = construct_transition_matrix(state_matrix_diseased_pre_prune, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_disease_pre_prune = pd.DataFrame(state_matrix_diseased_pre_prune)
    df_state_disease_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_disease_pre_prune.csv')
    df_transition_disease_pre_prune = pd.DataFrame(transition_matrix_diseased_pre_prune)
    df_transition_disease_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_disease_pre_prune.csv')
    #predict age
    predicted_ages_disease_pre_prune, actual_ages_disease_pre_prune = list_of_predicted_actual_age(df_disease,state_matrix_diseased_pre_prune,number_of_age_groups,list_of_methylation_sites,age_buckets_diseased_pre_prune)
    print('Diseased_Pre_Prune_Error', average_distance_apart(predicted_ages_disease_pre_prune, actual_ages_disease_pre_prune))
    graph_ages(predicted_ages_disease_pre_prune, actual_ages_disease_pre_prune,True)

    #Convert to Markov Model
    markov_chain_healthy_pre_prune = markov_chain(state_matrix_healthy_pre_prune, transition_matrix_healthy_pre_prune)
    df_markov_chain_healthy_pre_prune = pd.DataFrame(markov_chain_healthy_pre_prune)
    df_markov_chain_healthy_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_healthy_pre_prune.csv')
    markov_chain_disease_pre_prune = markov_chain(state_matrix_diseased_pre_prune, transition_matrix_diseased_pre_prune)
    df_markov_chain_disease_pre_prune = pd.DataFrame(markov_chain_disease_pre_prune)
    df_markov_chain_disease_pre_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_disease_pre_prune.csv')


    #Predict Healthy Age Post-Prune
    print('\n Post-Prune-Healthy_1')
    #make markov model
    state_matrix_healthy, age_buckets_healthy, methyl_buckets_healthy = construct_state_matrix_prune(df_healthy,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_healthy = construct_transition_matrix(state_matrix_healthy, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_healthy_post_prune = pd.DataFrame(state_matrix_healthy)
    df_state_healthy_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_healthy_post_prune.csv')
    df_transition_healthy_post_prune = pd.DataFrame(transition_matrix_healthy)
    df_transition_healthy_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_healthy_post_prune.csv')
    #predict age
    predicted_ages_healthy, actual_ages_healthy = list_of_predicted_actual_age(df_healthy,state_matrix_healthy,number_of_age_groups,list_of_methylation_sites,age_buckets_healthy)
    print('Healthy_Post_Prune_Error', average_distance_apart(predicted_ages_healthy, actual_ages_healthy))
    graph_ages(predicted_ages_healthy, actual_ages_healthy,False)

    #Predict Diseased Age Post-Prune 
    print('\n Post-Prune-Disease_1')
    #make markov model
    state_matrix_diseased, age_buckets_diseased, methyl_buckets_diseased = construct_state_matrix_prune(df_disease,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_diseased = construct_transition_matrix(state_matrix_diseased, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_disease_post_prune = pd.DataFrame(state_matrix_diseased)
    df_state_disease_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_disease_post_prune.csv')
    df_transition_disease_post_prune = pd.DataFrame(transition_matrix_diseased)
    df_transition_disease_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_disease_post_prune.csv')
    #predict age
    predicted_ages_disease, actual_ages_disease = list_of_predicted_actual_age(df_disease,state_matrix_diseased,number_of_age_groups,list_of_methylation_sites,age_buckets_diseased)
    print('Healthy_Pre_Prune_Error', average_distance_apart(predicted_ages_disease, actual_ages_disease))
    graph_ages(predicted_ages_disease, actual_ages_disease,True)

    #Convert Each State and Transition Matrix to a Markov Chain
    markov_chain_healthy_post_prune = markov_chain(state_matrix_healthy, transition_matrix_healthy)
    df_markov_chain_healthy_post_prune = pd.DataFrame(markov_chain_healthy_post_prune)
    df_markov_chain_healthy_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_healthy_post_prune.csv')
    markov_chain_disease_post_prune = markov_chain(state_matrix_diseased, transition_matrix_diseased)
    df_markov_chain_disease_post_prune = pd.DataFrame(markov_chain_disease_post_prune)
    df_markov_chain_disease_post_prune.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_disease_post_prune.csv')
    
    #condition Prediction Pre-Prune
    print('\n Condition Prediction')
    df_longitudinal_methylation = pd.read_csv("/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Input_Datasets/Longitude_Set.csv")
    predicted_condition, actual_condition = predictions_of_conditions(df_longitudinal_methylation, state_matrix_healthy_pre_prune, transition_matrix_healthy_pre_prune, state_matrix_diseased_pre_prune, transition_matrix_diseased_pre_prune, age_buckets_healthy_pre_prune, methyl_buckets_healthy, age_buckets_diseased_pre_prune, methyl_buckets_diseased)
    graph_conditions(predicted_condition,actual_condition)

    #Condition Prediction
    print('\n Condition Prediction')
    df_longitudinal_methylation = pd.read_csv("/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Input_Datasets/Longitude_Set.csv")
    predicted_condition, actual_condition = predictions_of_conditions(df_longitudinal_methylation, state_matrix_healthy, transition_matrix_healthy, state_matrix_diseased, transition_matrix_diseased, age_buckets_healthy, methyl_buckets_healthy, age_buckets_diseased, methyl_buckets_diseased)
    graph_conditions(predicted_condition,actual_condition)
    

    #Predict Healthy Age Post-Prune
    print('\n Post-Prune-Healthy_2')
    #make markov model
    state_matrix_healthy_2, age_buckets_healthy_2, methyl_buckets_healthy_2 = construct_state_matrix_prune_2(df_healthy,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_healthy_2 = construct_transition_matrix(state_matrix_healthy_2, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_healthy_post_prune_2 = pd.DataFrame(state_matrix_healthy_2)
    df_state_healthy_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_healthy_post_prune_2.csv')
    df_transition_healthy_post_prune_2 = pd.DataFrame(transition_matrix_healthy_2)
    df_transition_healthy_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_healthy_post_prune_2.csv')
    #predict age
    predicted_ages_healthy_2, actual_ages_healthy_2 = list_of_predicted_actual_age(df_healthy,state_matrix_healthy_2,number_of_age_groups,list_of_methylation_sites,age_buckets_healthy_2)
    print('Healthy_Post_Prune_Error', average_distance_apart(predicted_ages_healthy_2, actual_ages_healthy_2))
    graph_ages(predicted_ages_healthy_2, actual_ages_healthy_2,False)

    #Predict Diseased Age Post-Prune 
    print('\n Post-Prune-Disease_2')
    #make markov model
    state_matrix_diseased_2, age_buckets_diseased_2, methyl_buckets_diseased_2 = construct_state_matrix_prune_2(df_disease,list_of_methylation_sites, number_of_age_groups,number_of_methylation_sites,number_of_locations_at_each_age_group)
    transition_matrix_diseased_2 = construct_transition_matrix(state_matrix_diseased_2, number_of_age_groups, number_of_methylation_sites, number_of_locations_at_each_age_group)
    #save model to .csv file
    df_state_disease_post_prune_2 = pd.DataFrame(state_matrix_diseased_2)
    df_state_disease_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_state_disease_post_prune_2.csv')
    df_transition_disease_post_prune_2 = pd.DataFrame(transition_matrix_diseased_2)
    df_transition_disease_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_transition_disease_post_prune_2.csv')
    #predict age
    predicted_ages_disease_2, actual_ages_disease_2 = list_of_predicted_actual_age(df_disease,state_matrix_diseased_2,number_of_age_groups,list_of_methylation_sites,age_buckets_diseased_2)
    print('Healthy_Pre_Prune_Error', average_distance_apart(predicted_ages_disease_2, actual_ages_disease_2))
    graph_ages(predicted_ages_disease_2, actual_ages_disease_2,True)

    #Convert Each State and Transition Matrix to a Markov Chain
    markov_chain_healthy_post_prune_2 = markov_chain(state_matrix_healthy_2, transition_matrix_healthy_2)
    df_markov_chain_healthy_post_prune_2 = pd.DataFrame(markov_chain_healthy_post_prune_2)
    df_markov_chain_healthy_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_healthy_post_prune_2.csv')
    markov_chain_disease_post_prune_2 = markov_chain(state_matrix_diseased_2, transition_matrix_diseased_2)
    df_markov_chain_disease_post_prune_2 = pd.DataFrame(markov_chain_disease_post_prune_2)
    df_markov_chain_disease_post_prune_2.to_csv('/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Markov_Models/df_markov_chain_disease_post_prune_2.csv')
    
    #Condition Prediction
    print('\n Condition Prediction')
    df_longitudinal_methylation = pd.read_csv("/Users/raehash/Desktop/CMU Sophomore Year/Fall Semester/02-512 Computational Methods for Biological Modeling and Simulation/Markov_Model_Analysis/Input_Datasets/Longitude_Set.csv")
    predicted_condition_2, actual_condition_2 = predictions_of_conditions(df_longitudinal_methylation, state_matrix_healthy_2, transition_matrix_healthy_2, state_matrix_diseased_2, transition_matrix_diseased_2, age_buckets_healthy, methyl_buckets_healthy, age_buckets_diseased_2, methyl_buckets_diseased_2)
    graph_conditions(predicted_condition_2,actual_condition_2)


main()