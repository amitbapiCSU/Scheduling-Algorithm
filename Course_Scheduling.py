import os
import sys
from scipy.io import arff
import pandas
import codecs
import copy
import numpy as np
import random as rn
import matplotlib.pyplot as plt
import warnings

#./#Read_inputs.sh# instance1.crs instance1.stu instance1.sol

"""
tried code snippets
"""
"""
#f = open(sys.argv[1],"r")
#ft = f.read().decode("utf-8-sig").encode("utf-8")
#d = np.load(g)
#headers = f.readline().strip().split(DELIM)
#dt = np.load(f)
#d = np.loadtxt(ft, delimiter=DELIM)  # usecols=range(0,1)  ,delimiter='    '
#d = pandas.read_csv(open(sys.argv[1]))  # iris.data


#print d.values

#room_capacity = int(headers[0].strip())
#max_time_slots = int(headers[1])

for indx in range(no_of_courses) :
    course_conflict_matrix[indx, indx] = -1
        
course_list = student_course_map.get(1)
print int(course_list[0])

max_conflicted_courses = check_for_max_conflicted_courses(course_conflict_matrix)
time_slot = 1
day_slot = 1
crs1 = max_conflicted_courses[0]
crs2 = max_conflicted_courses[1]
timeslot_course_map[time_slot] = [crs1]
timeslot_course_map[time_slot+1] = [crs2]
day_course_map[day_slot] = [crs1, crs2]

day_1_exams = np.array(day_course_map.get(day_slot))
day_1_exams = day_1_exams.flatten()
day_course_map[day_slot] = day_1_exams

for timeslot in time_slots:
        courses_for_time_slot = time_slot_course_map.get(timeslot)
        for dayslot, courses in day_slot_course_map.items() :
            if np.array_equal(np.array(courses_for_time_slot), np.array(courses)) :
                for crs in courses :
                    if crs == 10 :
                        crs_str = '0'+str(crs)
                    else :
                        crs_str = '00'+str(crs)

                    updated_vals.append(crs_str+"    "+str(dayslot)+"    "+str(timeslot)+"\n")
                    break

"""

"""
Initialization of all variables
used in the assignment.
"""
DELIM = '    '
room_capacity = 0
max_time_slots = 0
every_day_time_slot = 5
tmp_arr = []
course_student_count_map = {}
student_course_map = {}
timeslot_course_map = {}
day_course_map = {}

course_time_slot_map = {}
course_day_slot_map = {}
neighbour_hood_soln = {}


temp_crs_map = {}
deluge_level = 0

"""
Configuration of conflict matrix which define for each pair
of couses #students common.
"""
def conflict_matrix_configuration(course_conflict_matrix, student_course_map, no_of_courses) :
    for indx in range(no_of_courses) :
            course_conflict_matrix[indx, indx] = -1
    course_conflict_matrix[:,0] = -2
        
    stu_ids = student_course_map.keys()
    for stu_id in stu_ids:
        courses = student_course_map.get(stu_id)
        #print range(len(courses))
        for indx1 in range(len(courses)) :
            for indx2 in range(len(courses)) :
                if int(indx2) > int (indx1) :
                    crs1 = int(courses[indx1])
                    crs2 = int(courses[indx2])
                    course_conflict_matrix[crs1, crs2] += 1
                    course_conflict_matrix[crs2, crs1] += 1
            #if indx == (len(courses) -1):
                #break
            #else :
            
    print course_conflict_matrix    
    return course_conflict_matrix
    
"""
Configuration of Initial Solution
"""
def config_initial_soln(course_conflict_matrix, course_student_count_map) :
    time_slot = 1
    day_slot = 1
    #day_course_map[day_slot] = []
    #day_course_map[day_slot+1] = []
    while True :
        #print 'sum '+str(np.sum(course_conflict_matrix[1:,0],axis = 0))
        if np.sum(course_conflict_matrix[1:,0],axis = 0) == -30 :
            break 
        max_conflicted_courses = check_for_max_conflicted_courses(course_conflict_matrix)
        crs1 = max_conflicted_courses[0]
        crs2 = max_conflicted_courses[1]
        #print 'conflicted crs '+str(crs1)+'   '+str(crs2)
        if course_conflict_matrix[crs1,0] == -2 :
            course_conflict_matrix[crs1,0] = -3
            timeslot_course_map[time_slot] = [crs1]    
            add_non_conflicted_courses(timeslot_course_map, time_slot, crs1, course_conflict_matrix)
            time_slot = time_slot+1
        
        if course_conflict_matrix[crs2,0] == -2 :
            course_conflict_matrix[crs2,0] = -3
            timeslot_course_map[time_slot] = [crs2]
            add_non_conflicted_courses(timeslot_course_map, time_slot, crs2, course_conflict_matrix)
            time_slot = time_slot+1
    
    timeslots = timeslot_course_map.keys()
    for timeslot in timeslots :
        #if timeslot <= 5:
        slot_key = int(len(timeslots)/every_day_time_slot)+1
        if day_course_map.has_key(slot_key) :
            day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
        else :
            day_course_map[slot_key] = []
            day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
        #else :
            #day_course_map[day_slot+1].append(timeslot_course_map.get(timeslot))
    
    print day_course_map
    return timeslot_course_map,day_course_map

"""
Determination of non-conflicting courses
"""
def add_non_conflicted_courses(timeslot_course_map, time_slot, crs, course_conflict_matrix) :
    crs_indx = 0
    if int(crs) == 10 :
        crs_key =  '0'+str(crs)
    else :
        crs_key = '00'+str(crs)
    print 'crs '+crs_key
    rm_capacity = course_student_count_map.get(crs_key)
    #print 'prev '+str(rm_capacity)
    for item in course_conflict_matrix[crs,:] :
        #print 'crs_indx '+str(crs_indx)
        #print 'item '+str(item)
        if item == 0 :
            if course_conflict_matrix[crs_indx,0] == -2 :
                if int(crs_indx) == 10 :
                    key =  '0'+str(crs_indx)
                else :
                    key = '00'+str(crs_indx)
                non_conflict_course_stud_cnt = course_student_count_map.get(key)
                #print 'non conflict crs '+key
                #print 'non conflict '+str(non_conflict_course_stud_cnt)
                prev_rm_cpcty = int(rm_capacity)
                rm_capacity = int(rm_capacity) + int(non_conflict_course_stud_cnt)
                #print 'combine rm '+str(rm_capacity)
                #print 'max '+str(room_capacity)
                if int(rm_capacity) <= int(room_capacity) :
                    timeslot_course_map[time_slot].append(crs_indx)
                    course_conflict_matrix[crs_indx,0] = -3
                    print timeslot_course_map
                else :
                    rm_capacity = prev_rm_cpcty
                #print course_conflict_matrix
        crs_indx = crs_indx+1
    
"""
Check for maximum conflicting courses
"""    
def check_for_max_conflicted_courses(course_conflict_matrix) :
    crs_indx = 1
    maxim = -6
    max_indices = ()
    while crs_indx< course_conflict_matrix.shape[0] :
        if course_conflict_matrix[crs_indx,0] == -2 :
            maxm = course_conflict_matrix[crs_indx,:].max()
            #print 'max '+str(maxm)            
            if maxim <= maxm :
                maxim = maxm
                max_col_indices = np.where(course_conflict_matrix[crs_indx,:] == course_conflict_matrix[crs_indx,:].max())
                max_indices = (crs_indx,max_col_indices[0][0])
        crs_indx =crs_indx+1
    return max_indices

"""
Total Cost Objective finction
"""
def student_cost_objective_func(student_course_map, timeslot_course_map, day_course_map) :
    timeslots = timeslot_course_map.keys()
    for timeslot in timeslots:
        courses = timeslot_course_map.get(timeslot)
        for course in courses:
            course_time_slot_map[course] = timeslot
        
    

    dayslots = day_course_map.keys()
    for dayslot in dayslots:
        course_arrs = day_course_map.get(dayslot)
        for course_arr in course_arrs :
            for indx in range(len(course_arr)):
                course_day_slot_map[course_arr[indx]] = dayslot

    
    students = student_course_map.keys()
    total_consecutive_penalty_arr = []
    total_overnight_penalty_arr = []
    
    for student in students :
        course_exams = student_course_map.get(student)
        #print  len(course_exams)
        overnightassignmentpenalty_stud = []
        consecutiveassignmentpenalty_stud = []

        for c1 in range(len(course_exams)):
            for c2 in range(len(course_exams)):
                if c1 != c2 :
                    day_slot_of_c1 =  course_day_slot_map.get(int(course_exams[c1]))
                    day_slot_of_c2 =  course_day_slot_map.get(int(course_exams[c2]))
        
                    time_slot_of_c1 = course_time_slot_map.get(int(course_exams[c1]))
                    time_slot_of_c2 = course_time_slot_map.get(int(course_exams[c2]))
        
                    if day_slot_of_c1 == day_slot_of_c2 : #consecutive
                        consecpenalty = 2 ** -abs(time_slot_of_c1 - time_slot_of_c2)
                        consecutiveassignmentpenalty_stud.append(consecpenalty)
                    else : #overnight
                        overnpenalty = 2 ** -abs(day_slot_of_c1 - day_slot_of_c2)
                        overnightassignmentpenalty_stud.append(overnpenalty)
        
        consecutiveassignmentpenaltysum_stud = np.sum(consecutiveassignmentpenalty_stud, axis = 0)
        overnightassignmentpenaltysum_stud = np.sum(overnightassignmentpenalty_stud, axis = 0)
        
        total_consecutive_penalty_arr.append(10*consecutiveassignmentpenaltysum_stud)
        total_overnight_penalty_arr.append(overnightassignmentpenaltysum_stud)
    
    total_consecutive_penalty = np.sum(total_consecutive_penalty_arr, axis = 0)
    total_overnight_penalty = np.sum(total_overnight_penalty_arr,axis = 0)
    
    total_student_cost = total_consecutive_penalty+total_overnight_penalty

    return total_student_cost

"""
Timeslot objective function
"""        
def time_slot_objective_function(timeslot_course_map) :
    timeslots = timeslot_course_map.keys()
    return len(timeslots)
    

def write_to_solution_file(file_name, course_time_slot_map, course_day_slot_map, cost_obj_val) :
    print 'in soln file'
    print course_time_slot_map
    updated_vals = []
    file_path = os.path.join('/s/chopin/k/grad/amchakra/Assignment2', file_name)
    time_slots = np.unique(course_time_slot_map.values())
    no_of_time_slots = len(time_slots)

    updated_vals.append(str(no_of_time_slots)+"    "+str(cost_obj_val)+"\n")
    
    courses = course_time_slot_map.keys()
    
    for crs in courses:
        if crs == 10 :
            crs_str = '0'+str(crs)
        else :
            crs_str = '00'+str(crs)
        
        timeslot = course_time_slot_map.get(crs)
        div = timeslot/every_day_time_slot
        if timeslot%every_day_time_slot == 0 :
            timeslot = timeslot - (div-1)*every_day_time_slot
        else :
            timeslot = timeslot - (div)*every_day_time_slot
        dayslot = course_day_slot_map.get(crs)
            
        updated_vals.append(crs_str+"    "+str(dayslot)+"    "+str(timeslot)+"\n")
    
    with open(file_path, "w") as fp:
        for indx in range(len(updated_vals)) :
            fp.write(str(updated_vals[indx]))
        #fp.write(updated_vals)

def get_neighbourhood_solution_space(timeslot_course_map) :
    neighbr_indx = 0
    for  key1 in timeslot_course_map.keys() :
        for key2 in timeslot_course_map.keys():
                
            if int(key2) > int(key1):
            
                key1_crs = []
                key2_crs = []
                temp_crs_map = copy.deepcopy(timeslot_course_map)
                        
            
                temp_crs_map[key1] = timeslot_course_map.get(key2)
                temp_crs_map[key2] = timeslot_course_map.get(key1)
            
                neighbour_hood_soln[neighbr_indx] = temp_crs_map
                neighbr_indx = neighbr_indx+1
        
    return neighbour_hood_soln

"""
The following code is to read the files passed
from comman line and parse them to poplulate 
corresponding variable.
"""

"""
loading *.crs file
"""
crs_file = codecs.open(sys.argv[1], "r", "utf-8")

itrn = 0
for line in crs_file:
    arr = line.split(DELIM)
    if itrn == 0 :
        room_capacity = arr[0]
        max_time_slots = arr[1]
        #print arr[0], arr[1]
    else :
        course_student_count_map[arr[0]] = arr[1]
    itrn = itrn+1


"""
loading *.stu file
"""

stu_file = codecs.open(sys.argv[2], "r", "utf-8")
student_id = 1
for line in stu_file:
    arr = line.split(DELIM)
    student_course_map[student_id] = arr
    student_id = student_id+1
#print sys.argv[2]

"""
Main section calling all functions
"""
no_of_courses = len(course_student_count_map.keys())+1
course_conflict_matrix = np.zeros((no_of_courses,no_of_courses))

course_conflict_matrix = conflict_matrix_configuration(course_conflict_matrix, student_course_map, no_of_courses)
 
indices = np.where(course_conflict_matrix[1,:] == course_conflict_matrix[1,:].max())

timeslot_course_map,day_course_map = config_initial_soln(course_conflict_matrix, course_student_count_map)

student_cost = student_cost_objective_func(student_course_map, timeslot_course_map, day_course_map)
no_of_timeslots = time_slot_objective_function(timeslot_course_map)

"""
Writing initial solution to a file
"""
write_to_solution_file(sys.argv[3], course_time_slot_map, course_day_slot_map, student_cost)

deluge_level = student_cost
deluge_level_decay_rate = 0.01 # need to be  passed as an input parameter
final_result_cost_func = int(student_cost/3)
no_iterations = int((student_cost - final_result_cost_func)/deluge_level_decay_rate)
print no_iterations

tmp_timeslot_course_map = copy.deepcopy(timeslot_course_map)
tmp_dayslot_course_map = {}
print timeslot_course_map
random_time_slot_schedule_selected = {}
day_slots_random_selected = {}
prev_student_cost = student_cost
write_file_count = 0

for indx in range(no_iterations) :
    print str(indx)+"  "+str(no_iterations)
    neighbour_hood_soln = get_neighbourhood_solution_space(tmp_timeslot_course_map)
    different_neighbour_hood_solns = neighbour_hood_soln.keys()
    rndom_index = rn.choice(different_neighbour_hood_solns)
    random_time_slot_schedule_selected = copy.deepcopy(neighbour_hood_soln.get(rndom_index))
    print random_time_slot_schedule_selected
    
    #day_slot = 1
    #day_slots_random_selected[day_slot] = []
    #day_slots_random_selected[day_slot+1] = []
    timeslots = random_time_slot_schedule_selected.keys()
    for timeslot in timeslots :
        #if timeslot <= 5:
        slot_key = int(len(timeslots)/every_day_time_slot)+1
        if day_slots_random_selected.has_key(slot_key) :
            day_slots_random_selected[slot_key].append(random_time_slot_schedule_selected.get(timeslot))
        else :
            day_slots_random_selected[slot_key] = []
            day_slots_random_selected[slot_key].append(random_time_slot_schedule_selected.get(timeslot))
    #for timeslot in timeslots :
        #if timeslot <= 5:
            #day_slots_random_selected[day_slot].append(random_time_slot_schedule_selected.get(timeslot))
        #else :
            #day_slots_random_selected[day_slot+1].append(random_time_slot_schedule_selected.get(timeslot))
    
    recent_student_cost = student_cost_objective_func(student_course_map, random_time_slot_schedule_selected, day_slots_random_selected)    
    
    if     recent_student_cost <= prev_student_cost or recent_student_cost <= deluge_level:
        write_file_count = 0
        prev_student_cost = recent_student_cost
        tmp_timeslot_course_map = copy.deepcopy(random_time_slot_schedule_selected)
        tmp_dayslot_course_map = copy.deepcopy(day_slots_random_selected)
    else :
        write_file_count = write_file_count+1
    
    if write_file_count == 5 :
        new_optim_student_cost = student_cost_objective_func(student_course_map, random_time_slot_schedule_selected, day_slots_random_selected)
        write_to_solution_file(sys.argv[3], course_time_slot_map, course_day_slot_map, new_optim_student_cost)
    
    deluge_level = deluge_level - deluge_level_decay_rate
    print 'prev_cost '+str(prev_student_cost)
    print 'new cost '+str(recent_student_cost) 


write_to_solution_file(sys.argv[3], course_time_slot_map, course_day_slot_map, prev_student_cost)

