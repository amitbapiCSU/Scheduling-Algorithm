import os
import sys
import codecs
import time
import copy
import math
import numpy as np
import random as rn

"""
Initialization of all variables
used in the assignment.
"""
DELIM = '	'
room_capacity = 0
max_time_slots = 0
every_day_time_slot = 5
tmp_arr = []
string_to_course_no_map = {}
course_no_to_string_map = {}
course_student_count_map = {}
student_course_map = {}
timeslot_course_map = {}
day_course_map = {}

course_time_slot_map = {}
course_day_slot_map = {}
neighbour_hood_soln = {}

initial_population_indices = {}
initial_course_timeslot_population = {}
initial_course_dayslot_population = {}
initial_timeslot_course_population = {}
initial_dayslot_course_population = {}

next_population_indices = {}
next_course_timeslot_population = {}
next_course_dayslot_population = {}
next_timeslot_course_population = {}
next_dayslot_course_population = {}

var_population_indices = {}
var_course_timeslot_population = {}
var_course_dayslot_population = {}
var_timeslot_course_population = {}
var_dayslot_course_population = {}

fitness_scores_by_population = {}
fitness_scores_by_population_div = {}

combine_prob = 0.7
mutation_prob = 0.003

temp_crs_map = {}
tmp_day_course_map = {}
tmp_course_time_slot_map = {}
tmp_course_day_slot_map = {}

deluge_level = 0

student_cost = 0
no_of_timeslots = 0
min_no_of_timeslots = 0

optim_timeslot_or_cost = sys.argv[4]

no_of_crs = 50
min_req_crs = 2
max_req_crs = 6
no_of_studs = 50
stu_crs_map = {}
courses_flag = {}
crs_stud_cnt = {}

file_cnt = 0
write_file_count = 0

"""
Configuration of conflicting matrix which define for each pair
of courses the number of students common.
"""
def conflict_matrix_configuration(course_conflict_matrix, student_course_map, no_of_courses) :
	for indx in range(no_of_courses) :
			course_conflict_matrix[indx, indx] = -1
	course_conflict_matrix[:,0] = -2 # the first column is made -2
	
	# computation of number of common students for each pair of courses	
	stu_ids = student_course_map.keys()
	for stu_id in stu_ids:
		courses = student_course_map.get(stu_id)
		for indx1 in range(len(courses)) :
			for indx2 in range(len(courses)) :
				if int(indx2) > int (indx1) :
					crs1 = string_to_course_no_map.get(courses[indx1])
					#print courses
					#print indx2
					crs2 = string_to_course_no_map.get(courses[indx2])
					course_conflict_matrix[crs1, crs2] += 1
					course_conflict_matrix[crs2, crs1] += 1
	#print course_conflict_matrix		
	return course_conflict_matrix
	
"""
Configuration of Initial Solution, which is meant to provide
a schedule of minimum number of timeslots.
"""
def config_initial_soln(course_conflict_matrix, course_student_count_map) :
	
	time_slot = 1
	day_slot = 1

	while True :
		if np.sum(course_conflict_matrix[1:,0],axis = 0) == -3*int(no_of_courses-1) :
			break 
		max_conflicted_courses = check_for_max_conflicted_courses(course_conflict_matrix)
		#print max_conflicted_courses
		crs1 = max_conflicted_courses[0]
		crs2 = max_conflicted_courses[1]
		
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
		slot_key = 0
		if timeslot%every_day_time_slot == 0 :
			slot_key = int(timeslot/every_day_time_slot)
		else :
			slot_key = int(timeslot/every_day_time_slot)+1
			
		if day_course_map.has_key(slot_key) :
			day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
		else :
			day_course_map[slot_key] = []
			day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
	
	for timeslot in timeslots :
		crss = timeslot_course_map.get(timeslot)
		cap = 0
		#print "courses of slot "+str(timeslot)
		#print crss
		for crs in crss :
			cap = int(cap) + int(course_student_count_map.get(crs))
		#print "slot "+str(timeslot)+" capacity "+str(cap)
		
	return timeslot_course_map,day_course_map

"""
Construction of initial solution and computation 
of initial student costobjective function
"""
def construct_initial_solution_compute_obj_val(course_conflict_matrix, student_course_map, no_of_courses, course_student_count_map, timeslot_course_map, day_course_map) :
	course_conflict_matrix = conflict_matrix_configuration(course_conflict_matrix, student_course_map, no_of_courses)
	timeslot_course_map,day_course_map = config_initial_soln(course_conflict_matrix, course_student_count_map)
	course_time_slot_map,course_day_slot_map = construct_course_timeslot_and_dayslot_map(timeslot_course_map, day_course_map)
	student_cost = student_cost_objective_func(student_course_map, course_time_slot_map, course_day_slot_map)
	return student_cost

"""
This function checks conflicts of the currently selected course 
with the scheduled courses in that timeslot 
"""
def check_conflict_with_existing_timeslot_courses(crs_indx, course_conflict_matrix, time_slot) :
	existing_scheduled_courses = timeslot_course_map.get(time_slot)
	cnfct_flg = False
	for exst_crs in existing_scheduled_courses:
		 if course_conflict_matrix[crs_indx,exst_crs] != 0 :
			 cnfct_flg = True
			 break
	return cnfct_flg
			 
"""
Determination of non-conflicting courses
"""
def add_non_conflicted_courses(timeslot_course_map, time_slot, crs, course_conflict_matrix) :
	crs_indx = 0
	crskey = crs #crs_key(crs)
	rm_capacity = course_student_count_map.get(crskey)
	
	for item in course_conflict_matrix[crs,:] :
		cnflict_flag = False
		items = timeslot_course_map.get(time_slot)
		if item == 0 :
			
			if len(items) > 1 :
				cnflict_flag = check_conflict_with_existing_timeslot_courses(crs_indx, course_conflict_matrix, time_slot)
				
			if course_conflict_matrix[crs_indx,0] == -2 and cnflict_flag == False :
				key = crs_indx #crs_key(int(crs_indx))
				#print "key "+str(key)
				non_conflict_course_stud_cnt = course_student_count_map.get(key)
				prev_rm_cpcty = int(rm_capacity)
				tmp_rm_capacity = int(rm_capacity) + int(non_conflict_course_stud_cnt)
				
				if int(tmp_rm_capacity) <= int(room_capacity) :
					timeslot_course_map[time_slot].append(crs_indx)
					course_conflict_matrix[crs_indx,0] = -3
					rm_capacity = tmp_rm_capacity
				#else :
					#rm_capacity = prev_rm_cpcty
		crs_indx = crs_indx+1
	#print "crs "+str(crs)+" crs rm cpty "+str(rm_capacity)+" max cap "+str(room_capacity)
	
"""
Check for maximum conflicting courses
"""	
def check_for_max_conflicted_courses(course_conflict_matrix) :
	crs_indx = 1
	maxim = -6
	max_indices = ()
	while crs_indx < course_conflict_matrix.shape[0] :
		if course_conflict_matrix[crs_indx,0] == -2 :
			maxm = course_conflict_matrix[crs_indx,:].max()	
			if maxim <= maxm :
				maxim = maxm
				max_col_indices = np.where(course_conflict_matrix[crs_indx,:] == course_conflict_matrix[crs_indx,:].max())
				#print 'crs '+str(crs_indx)+' maxm '+str(maxm)+' col indices '+str(max_col_indices)
				max_indices = (crs_indx,max_col_indices[0][0])
		crs_indx =crs_indx+1
	#print 'max indices '+str(max_indices)
	return max_indices


"""
Computation of Total student Cost Objective function
"""	
def student_cost_objective_func(student_course_map, crs_timeslot_map, crs_day_map) :
	
	#course_time_slot_map, course_day_slot_map = construct_course_timeslot_and_dayslot_map(timeslot_course_map, day_course_map)
	#print "in objective"
	#print crs_timeslot_map
	#print crs_day_map
	#print string_to_course_no_map
	
	students = student_course_map.keys()
	total_consecutive_penalty_arr = []
	total_overnight_penalty_arr = []
	
	for student in students :
		crs_exams = student_course_map.get(student)
		for c1 in range(len(crs_exams)):
			for c2 in range(len(crs_exams)):
				if c1 != c2 :
					if crs_exams[c1] == crs_exams[c2] :
						print "same course multiple times for a single student "
						#print "c1 "+crs_exams[c1]
						#print "c2 "+crs_exams[c2]
						print "For student id "+str(student)+" same courses ("+crs_exams[c1]+","+crs_exams[c2]+")"
						break
					
	for student in students :
		course_exams = student_course_map.get(student)
		overnightassignmentpenalty_stud = []
		consecutiveassignmentpenalty_stud = []

		for c1 in range(len(course_exams)):
			for c2 in range(len(course_exams)):
				if c1 != c2 :
					#print course_exams[c1]
					#print str(string_to_course_no_map.get(course_exams[c1]))
					#print course_exams[c2]
					#print str(string_to_course_no_map.get(course_exams[c2]))
					day_slot_of_c1 =  crs_day_map.get(string_to_course_no_map.get(course_exams[c1]))
					day_slot_of_c2 =  crs_day_map.get(string_to_course_no_map.get(course_exams[c2]))
		
					time_slot_of_c1 = crs_timeslot_map.get(string_to_course_no_map.get(course_exams[c1]))
					#print "time 1 1 -> "+str(time_slot_of_c1)
					#print "day 1 1 -> "+str(day_slot_of_c1)
					#timeslot = course_time_slot_map.get(crs)
					
					div = time_slot_of_c1/every_day_time_slot
					if time_slot_of_c1%every_day_time_slot == 0 :
						time_slot_of_c1 = time_slot_of_c1 - (div-1)*every_day_time_slot
					else :
						time_slot_of_c1 = time_slot_of_c1 - (div)*every_day_time_slot
					#print "time 1 2 -> "+str(time_slot_of_c1)	
					
					time_slot_of_c2 = crs_timeslot_map.get(string_to_course_no_map.get(course_exams[c2]))
					#print "time 2 1 -> "+str(time_slot_of_c2)
					#print "day 2 1 -> "+str(day_slot_of_c2)
					
					div = time_slot_of_c2/every_day_time_slot
					if time_slot_of_c2%every_day_time_slot == 0 :
						time_slot_of_c2 = time_slot_of_c2 - (div-1)*every_day_time_slot
					else :
						time_slot_of_c2 = time_slot_of_c2 - (div)*every_day_time_slot
					#print "time 2 2 -> "+str(time_slot_of_c2)
					
					if day_slot_of_c1 == day_slot_of_c2 : #consecutive
						consecpenalty = 2 ** -abs(time_slot_of_c1 - time_slot_of_c2)
						consecutiveassignmentpenalty_stud.append(consecpenalty)
					else : #overnight
						overnpenalty = 2 ** -abs(day_slot_of_c1 - day_slot_of_c2)
						overnightassignmentpenalty_stud.append(overnpenalty)
		
		consecutiveassignmentpenaltysum_stud = np.sum(consecutiveassignmentpenalty_stud, axis = 0)
		overnightassignmentpenaltysum_stud = np.sum(overnightassignmentpenalty_stud, axis = 0)
		#print "courses stud "
		#print course_exams
		#print " consec "+str(consecutiveassignmentpenaltysum_stud)+" overn "+str(overnightassignmentpenaltysum_stud)
		total_consecutive_penalty_arr.append(10*consecutiveassignmentpenaltysum_stud)
		total_overnight_penalty_arr.append(overnightassignmentpenaltysum_stud)
	
	total_consecutive_penalty = np.sum(total_consecutive_penalty_arr, axis = 0)
	total_overnight_penalty = np.sum(total_overnight_penalty_arr,axis = 0)
	
	total_student_cost = total_consecutive_penalty+total_overnight_penalty

	return total_student_cost

"""
Construction of course timeslot map and course dayslot map
"""
def construct_course_timeslot_and_dayslot_map(timeslot_course_map, day_course_map) :
	#print timeslot_course_map
	#print "2"
	#print day_course_map
	
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
	
	return course_time_slot_map, course_day_slot_map
	
"""
This functions gets neighbourhood of a schedule by
swapping any two timeslots of the schedule
"""
def get_neighbourhood_solution_space(timeslot_course_map) :
	neighbr_indx = 0
	for  key1 in timeslot_course_map.keys() :
		for key2 in timeslot_course_map.keys():
				
			if int(key2) > int(key1):
				temp_crs_map = copy.deepcopy(timeslot_course_map)
						
				temp_crs_map[key1] = timeslot_course_map.get(key2)
				temp_crs_map[key2] = timeslot_course_map.get(key1)
			
				neighbour_hood_soln[neighbr_indx] = temp_crs_map
				neighbr_indx = neighbr_indx+1
		
	return neighbour_hood_soln


def get_all_relevant_maps(temp_crs_map, tmp_day_course_map,tmp_course_time_slot_map, tmp_course_day_slot_map) :
	timeslots = temp_crs_map.keys()	
	
	for timeslot in timeslots :
		slot_key = 0
		if timeslot%every_day_time_slot == 0 :
			slot_key = int(timeslot/every_day_time_slot)
		else :
			slot_key = int(timeslot/every_day_time_slot)+1
			
		if tmp_day_course_map.has_key(slot_key) :
			tmp_day_course_map[slot_key].append(temp_crs_map.get(timeslot))
		else :
			tmp_day_course_map[slot_key] = []
			tmp_day_course_map[slot_key].append(temp_crs_map.get(timeslot))
	
	tmp_course_time_slot_map, tmp_course_day_slot_map = construct_course_timeslot_and_dayslot_map(temp_crs_map, tmp_day_course_map)
	
	return temp_crs_map, tmp_day_course_map,tmp_course_time_slot_map, tmp_course_day_slot_map
"""
Configuration of initial population
"""	
def config_initial_population(timeslot_course_map, day_course_map) :
	population_indx = 0
	course_time_slot_map, course_day_slot_map = construct_course_timeslot_and_dayslot_map(timeslot_course_map, day_course_map)
	initial_population_indices[population_indx] = population_indx
	initial_course_timeslot_population[population_indx] = course_time_slot_map
	initial_course_dayslot_population[population_indx] = course_day_slot_map
	initial_timeslot_course_population[population_indx] = timeslot_course_map
	initial_dayslot_course_population[population_indx] = day_course_map
	
	#print "total loop "+str(2*(len(timeslot_course_map.keys())))
	cnt = 0
	for  key1 in timeslot_course_map.keys() :
		for key2 in timeslot_course_map.keys():
			#print "in config population "	+str(cnt)
			if (int(cnt) > 1600) : #maximum initial population size
				break
			cnt = cnt + 1
			if int(key2) > int(key1):
				temp_crs_map = copy.deepcopy(timeslot_course_map)
						
				temp_crs_map[key1] = timeslot_course_map.get(key2)
				temp_crs_map[key2] = timeslot_course_map.get(key1)
				
				tmp_day_course_map.clear()
				#if tmp_course_time_slot_map.haskeys() :
					#tmp_course_time_slot_map.clear()
					#tmp_course_day_slot_map.clear()
					
				timeslots = temp_crs_map.keys()	
	
				for timeslot in timeslots :
					slot_key = 0
					if timeslot%every_day_time_slot == 0 :
						slot_key = int(timeslot/every_day_time_slot)
					else :
						slot_key = int(timeslot/every_day_time_slot)+1
			
					if tmp_day_course_map.has_key(slot_key) :
						tmp_day_course_map[slot_key].append(temp_crs_map.get(timeslot))
					else :
						tmp_day_course_map[slot_key] = []
						tmp_day_course_map[slot_key].append(temp_crs_map.get(timeslot))
				
				#print "tmp"
				#print temp_crs_map
				#print tmp_day_course_map
				
				tmp_course_time_slot_map, tmp_course_day_slot_map = construct_course_timeslot_and_dayslot_map(temp_crs_map, tmp_day_course_map)
				
				#print "tmp time"
				#print tmp_course_time_slot_map
				#print tmp_course_day_slot_map
				population_indx = population_indx+1
				
				#temp_crs_map, tmp_day_course_map, tmp_course_time_slot_map, tmp_course_day_slot_map = get_all_relevant_maps(temp_crs_map, tmp_day_course_map, tmp_course_time_slot_map, tmp_course_day_slot_map) 
				
				initial_population_indices[population_indx] = population_indx
				initial_course_timeslot_population[population_indx] = copy.deepcopy(tmp_course_time_slot_map)
				initial_course_dayslot_population[population_indx] = copy.deepcopy(tmp_course_day_slot_map)
				initial_timeslot_course_population[population_indx] = copy.deepcopy(temp_crs_map)
				initial_dayslot_course_population[population_indx] = copy.deepcopy(tmp_day_course_map)	
							
				#print "pop"
				#print initial_course_timeslot_population
				#neighbour_hood_soln[neighbr_indx] = temp_crs_map
				#neighbr_indx = neighbr_indx+1
	
"""
This function moves timeslots in a schedule to make the schedule 
sparse and optimize student cost
"""			
def get_neighbourhood_solution_space_1(timeslot_course_map) :
	#print ('in diff neighbor')
	#print max_time_slots
	
	timeslots = timeslot_course_map.keys()
	print timeslots
	timeslot_course_map_len = len(timeslots)
	print timeslot_course_map_len
	tmp_day_course_map = {}
	
	while (int(timeslot_course_map_len) < int(max_time_slots)) :
		#print str(timeslot_course_map_len)+"	"+str(max_time_slots)
		timeslots = timeslot_course_map.keys()
		for timeslot in timeslots:
			course_for_timeslot  = timeslot_course_map.get(timeslot)
			max_stud_count_fr_timeslot = 0
			crs_wid_max_stud = 0
			total_stud_count_fr_timeslot = 0
			if len(course_for_timeslot) >= 2 :
				for crs in course_for_timeslot :
					crs_str = crs #crs_key(crs)
					crs_stud_cnt = course_student_count_map.get(crs_str)
				
					if crs_stud_cnt > max_stud_count_fr_timeslot:
						max_stud_count_fr_timeslot = crs_stud_cnt
						crs_wid_max_stud = crs
				
					total_stud_count_fr_timeslot = int(total_stud_count_fr_timeslot) + int(crs_stud_cnt)
				
				if int(total_stud_count_fr_timeslot) <= int(room_capacity) and int(timeslot_course_map_len) < int(max_time_slots) :
					#print "in spread "+str(timeslot_course_map_len)
					course_for_timeslot.remove(crs_wid_max_stud)
					timeslot_course_map[timeslot] = [crs_wid_max_stud]
					timeslot_course_map[timeslot_course_map_len+1] = course_for_timeslot
					timeslot_course_map_len = len(timeslot_course_map.keys())
					
		#break
	#print "assigning day slot"
	day_course_map.clear()
	timeslots = timeslot_course_map.keys()	
	for timeslot in timeslots :
		slot_key = 0
		if timeslot%every_day_time_slot == 0 :
			slot_key = int(timeslot/every_day_time_slot)
		else :
			slot_key = int(timeslot/every_day_time_slot)+1
		#print "timeslot "+str(timeslot)+"	dayslot "+str(slot_key)	
		if day_course_map.has_key(slot_key) :
			day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
		else :
			day_course_map[slot_key] = []
			day_course_map[slot_key].append(timeslot_course_map.get(timeslot))
	
	print timeslot_course_map	
	return timeslot_course_map,day_course_map
	#return get_neighbourhood_solution_space(timeslot_course_map)
					
"""
This function writes the solution (.sol) file
"""
def write_to_solution_file(file_name, course_time_slot_map, course_day_slot_map, cost_obj_val) :
	updated_vals = []
	file_path = os.path.join(os.getcwd(), file_name)
	time_slots = np.unique(course_time_slot_map.values())
	no_of_time_slots = len(time_slots)

	updated_vals.append(str(no_of_time_slots)+"	"+str(cost_obj_val)+"\n")
	
	courses = course_time_slot_map.keys()
	
	for crs in courses:
		crs_str = course_no_to_string_map.get(crs)#crs_key(crs)
		
		timeslot = course_time_slot_map.get(crs)
		div = timeslot/every_day_time_slot
		if timeslot%every_day_time_slot == 0 :
			timeslot = timeslot - (div-1)*every_day_time_slot
		else :
			timeslot = timeslot - (div)*every_day_time_slot
			
		dayslot = course_day_slot_map.get(crs)
			
		updated_vals.append(str(crs_str)+"	"+str(dayslot)+"    "+str(timeslot)+"\n")
		
	with open(file_path, "w") as fp:
		for indx in range(len(updated_vals)) :
			fp.write(str(updated_vals[indx]))	
"""
Fitness Score Calculation
"""
def fitness_score_calc(var_population_indices, var_course_timeslot_population, var_course_dayslot_population) :
	populations = var_population_indices.keys()
	total_fitness = 0
	print "length of current population "+str(len(populations))
	for elem in populations :
		print "elem "+str(elem)
		stu_cost = student_cost_objective_func(student_course_map, var_course_timeslot_population.get(elem), var_course_dayslot_population.get(elem))
		fitness_scores_by_population[elem] = 1/(1+stu_cost)**2
		total_fitness = total_fitness + fitness_scores_by_population.get(elem)
	
	#for elem in populations :
		#fitness_scores_by_population_div[elem] = fitness_scores_by_population.get(elem)/total_fitness
	
	return total_fitness

"""
Roulette Wheel selection
"""
def roulette_wheel_selection(var_population_indices, total_fitness) :
	populations = var_population_indices.keys()
	
	randm_no = rn.uniform(0,total_fitness)  # range(int(total_fitness)+1)
	#print "rndm no "+str(randm_no)
	summ = 0
	iselected = 2
	#while(True) :
	for elem in populations :
		#print "summ "+str(summ)
		if (summ > randm_no):
			#iselected = elem - 1
			break
		summ = summ + fitness_scores_by_population.get(elem)
		iselected = elem
	#print "selected "+str(iselected)
	return iselected
	
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
	arr = line.split()
	#print arr
	if itrn == 0 :
		room_capacity = arr[0]
		max_time_slots = arr[1]
		#print arr[0], arr[1]
	else :
		course_no_to_string_map[itrn] = arr[0]
		string_to_course_no_map[arr[0]] = itrn
		course_student_count_map[itrn] = arr[1]
	itrn = itrn+1


"""
loading *.stu file
"""
stu_file = codecs.open(sys.argv[2], "r", "utf-8")
student_id = 1
for line in stu_file:
	arr = line.split()
	student_course_map[student_id] = arr
	student_id = student_id+1
	
no_of_courses = len(course_student_count_map.keys())+1
course_conflict_matrix = np.zeros((no_of_courses,no_of_courses))

course_conflict_matrix = conflict_matrix_configuration(course_conflict_matrix, student_course_map, no_of_courses)
timeslot_course_map,day_course_map = config_initial_soln(course_conflict_matrix, course_student_count_map)
print timeslot_course_map
min_no_of_timeslots = int(len(timeslot_course_map.keys()))
if int(optim_timeslot_or_cost) == 1 :
	timeslot_course_map,day_course_map = get_neighbourhood_solution_space_1(timeslot_course_map)

config_initial_population(timeslot_course_map, day_course_map)
"""
Writing initial solution to file and
track the student cost to optimize
"""
print "to write initial sol"
corse_time_slot_map, corse_day_slot_map = construct_course_timeslot_and_dayslot_map(timeslot_course_map, day_course_map) 
prev_stu_cost = student_cost_objective_func(student_course_map, corse_time_slot_map, corse_day_slot_map)
write_to_solution_file(sys.argv[3], corse_time_slot_map, corse_day_slot_map, prev_stu_cost)

#t1 = time.clock()
print "writing done"

"""
Genteic Algorithm
"""

itrtn = 0
while(True) :  # it is expected to stop by 10 minutes externally by the checker #(time.clock() - t1) < 600
	print "init population"
	#print initial_population_indices
	if (itrtn == 0) :
		var_population_indices = copy.deepcopy(initial_population_indices)
		var_course_timeslot_population = copy.deepcopy(initial_course_timeslot_population)
		var_course_dayslot_population = copy.deepcopy(initial_course_dayslot_population)
		var_timeslot_course_population = copy.deepcopy(initial_timeslot_course_population)
		var_dayslot_course_population = copy.deepcopy(initial_dayslot_course_population)
	else  :
		var_population_indices = copy.deepcopy(next_population_indices)
		var_course_timeslot_population = copy.deepcopy(next_course_timeslot_population)
		var_course_dayslot_population = copy.deepcopy(next_course_dayslot_population)
		var_timeslot_course_population = copy.deepcopy(next_timeslot_course_population)
		var_dayslot_course_population = copy.deepcopy(next_dayslot_course_population)
		
		next_population_indices.clear()
		next_course_timeslot_population.clear()
		next_course_dayslot_population.clear()
		next_timeslot_course_population.clear()
		next_dayslot_course_population.clear()
	
	itrtn = itrtn + 1
	print "iterations "+str(itrtn)
	#print "var population"
	#print var_population_indices
	
	total_fitness = fitness_score_calc(var_population_indices, var_course_timeslot_population, var_course_dayslot_population)

	#total_fitness_1 = np.sum(np.array(fitness_scores_by_population.values()))

	#print "1"
	#print fitness_scores_by_population
	populations = var_population_indices.keys()
	maxim = 0.0
	imax = 0
	for elem in populations :
		if fitness_scores_by_population.get(elem) > maxim:
			maxim = fitness_scores_by_population.get(elem)
			imax = elem
	
	print "fittest "+str(imax)
	
	stu_cost = student_cost_objective_func(student_course_map, var_course_timeslot_population.get(imax), var_course_dayslot_population.get(imax))
	print "stu cost "+str(stu_cost)
	if int(stu_cost) <= int(prev_stu_cost) :
		write_to_solution_file(sys.argv[3], var_course_timeslot_population.get(imax), var_course_dayslot_population.get(imax), stu_cost)
		prev_stu_cost = stu_cost
	
	#To do get the fittest
	#check whether the solution is okay
	#check for initial sol if necessary
	
	population_indx = 0
	while(population_indx <= len(populations)) : 
		prnt_chrsm_1_indx = roulette_wheel_selection(var_population_indices, total_fitness)
		prnt_chrsm_2_indx = roulette_wheel_selection(var_population_indices, total_fitness)
		
		prnt_chromsm_1_crs_time = copy.deepcopy(var_course_timeslot_population.get(prnt_chrsm_1_indx))
		prnt_chromsm_2_crs_time = copy.deepcopy(var_course_timeslot_population.get(prnt_chrsm_2_indx))

		prnt_chromsm_1_crs_day = copy.deepcopy(var_course_dayslot_population.get(prnt_chrsm_1_indx))
		prnt_chromsm_2_crs_day = copy.deepcopy(var_course_dayslot_population.get(prnt_chrsm_2_indx))

		prnt_chromsm_1_time_crs = copy.deepcopy(var_timeslot_course_population.get(prnt_chrsm_1_indx))
		prnt_chromsm_2_time_crs = copy.deepcopy(var_timeslot_course_population.get(prnt_chrsm_2_indx))

		prnt_chromsm_1_day_crs = copy.deepcopy(var_dayslot_course_population.get(prnt_chrsm_1_indx))
		prnt_chromsm_2_day_crs = copy.deepcopy(var_dayslot_course_population.get(prnt_chrsm_2_indx))

		

		#print "parent day "
		#print prnt_chromsm_1_crs_day
		#print prnt_chromsm_2_crs_day
		#combine_prob = 0.7
		#mutation_prob = 0.003
		child_chrsm_1_crs_time = {}
		child_chrsm_1_time_crs = {}
		child_chrsm_2_crs_time = {}
		child_chrsm_2_time_crs = {}
		child_chrsm_1_crs_day = {}
		child_chrsm_1_day_crs= {}
		child_chrsm_2_crs_day = {}
		child_chrsm_2_day_crs = {}
		select_child_1 = True
		select_child_2 = True

		
		#Cross Over (fpu) on two parent chrosomes
		
		if rn.choice(range(0,1)) < combine_prob :
	
			fixed_nos = 2
			rndm_fixed_nos = rn.sample(range(1,len(prnt_chromsm_1_crs_time.keys())+1), fixed_nos)
	
			for crs in prnt_chromsm_1_crs_time.keys() :
				if crs in rndm_fixed_nos :
					
					child_chrsm_1_crs_time[crs] = prnt_chromsm_1_crs_time.get(crs)
					child_chrsm_1_crs_day[crs] = prnt_chromsm_1_crs_day.get(crs)
					child_chrsm_2_crs_time[crs] = prnt_chromsm_2_crs_time.get(crs)
					child_chrsm_2_crs_day[crs] = prnt_chromsm_2_crs_day.get(crs)
			
				else :
					child_chrsm_1_crs_time[crs] = prnt_chromsm_2_crs_time.get(crs)
					child_chrsm_1_crs_day[crs] = prnt_chromsm_2_crs_day.get(crs)
					child_chrsm_2_crs_time[crs] = prnt_chromsm_1_crs_time.get(crs)
					child_chrsm_2_crs_day[crs] = prnt_chromsm_1_crs_day.get(crs)
		else :
			child_chrsm_1_crs_time = copy.deepcopy(prnt_chromsm_1_crs_time)
			child_chrsm_1_crs_day = copy.deepcopy(prnt_chromsm_1_crs_day)
			child_chrsm_2_crs_time = copy.deepcopy(prnt_chromsm_2_crs_time)
			child_chrsm_2_crs_day = copy.deepcopy(prnt_chromsm_2_crs_day)

		
		#mutation on each child (random swappping of two courses with 
		#corresponding time and day slot)
		
		if rn.choice(range(0,1)) < mutation_prob :
			courses = child_chrsm_1_crs_time.keys()
			rndm_seleted_crss = rn.sample(courses, 2)
	
			crs1 = rndm_seleted_crss[0]
			crs2 = rndm_seleted_crss[1]
	
			tmp_time_slot = child_chrsm_1_crs_time.get(crs1)
			child_chrsm_1_crs_time[crs1] = child_chrsm_1_crs_time.get(crs2)
			child_chrsm_1_crs_time[crs2] = tmp_time_slot
	
			tmp_day_slot = child_chrsm_1_crs_day.get(crs1)
			child_chrsm_1_crs_day[crs1] = child_chrsm_1_crs_day.get(crs2)
			child_chrsm_1_crs_day[crs2] = tmp_day_slot
	
			courses = child_chrsm_2_crs_time.keys()
			rndm_seleted_crss = rn.sample(courses, 2)
	
			crs1 = rndm_seleted_crss[0]
			crs2 = rndm_seleted_crss[1]
	
			tmp_time_slot = child_chrsm_2_crs_time.get(crs1)
			child_chrsm_2_crs_time[crs1] = child_chrsm_2_crs_time.get(crs2)
			child_chrsm_2_crs_time[crs2] = tmp_time_slot
	
			tmp_day_slot = child_chrsm_2_crs_day.get(crs1)
			child_chrsm_2_crs_day[crs1] = child_chrsm_2_crs_day.get(crs2)
			child_chrsm_2_crs_day[crs2] = tmp_day_slot

		## discarding and decreasing left
		#print course_conflict_matrix
		
		#If infeasible child is produced 
		#after crsossover and mutation then
		#rather repairing just discard it. 
		#child 1
		
		for crs in 	child_chrsm_1_crs_time.keys() :
			slot = child_chrsm_1_crs_time.get(crs)
			if child_chrsm_1_time_crs.has_key(slot) :
				arr = child_chrsm_1_time_crs.get(slot)
				arr.append(crs)
				child_chrsm_1_time_crs[slot] = arr
			else :
				child_chrsm_1_time_crs[slot] = []
				child_chrsm_1_time_crs[slot] = [crs]
	
		timeslots = child_chrsm_1_time_crs.keys()	
		
		for slot in timeslots:
			crss = child_chrsm_1_time_crs.get(slot)
			#print "timeslots "+str(len(timeslots))+" max slots "+str(int(max_time_slots))+" min slots "+str(int(min_no_of_timeslots))
			if int(optim_timeslot_or_cost) == 1 and int(len(timeslots)) != int(max_time_slots) :
				select_child_1 = False
				break
			if int(optim_timeslot_or_cost) == 0 and int(len(timeslots)) != int(min_no_of_timeslots) :
				select_child_1 = False
				break
			if len(crss) > 1 :
				for indx1 in range(0,len(crss)) :
					for indx2 in range(0,len(crss)) :
						if int(indx1) != int(indx2) :
							crs1 = int(crss[indx1])
							crs2 = int(crss[indx2])
							#print "1 course conflict "+str(crs1)+"	"+str(crs2)+" val "+str(int(course_conflict_matrix[crs1,crs2]))
							if int(course_conflict_matrix[crs1,crs2]) > 0 :
								select_child_1 = False
								break
					#break
			#break
			if select_child_1 :
				capacity = 0
				for crs in crss :
					capacity = int(capacity) + int(course_student_count_map.get(crs))
		
				if int(capacity) > int(room_capacity):
					select_child_1 = False
					break
		
		
		#If child 1 produced as a feasible one then populate it to next population
		
		if select_child_1 :
			for timeslot in timeslots :
				slot_key = 0
				if timeslot%every_day_time_slot == 0 :
					slot_key = int(timeslot/every_day_time_slot)
				else :
					slot_key = int(timeslot/every_day_time_slot)+1
			
				if child_chrsm_1_day_crs.has_key(slot_key) :
					child_chrsm_1_day_crs[slot_key].append(child_chrsm_1_time_crs.get(timeslot))
				else :
					child_chrsm_1_day_crs[slot_key] = []
					child_chrsm_1_day_crs[slot_key].append(child_chrsm_1_time_crs.get(timeslot))
			
			next_population_indices[population_indx] = population_indx
			next_course_timeslot_population[population_indx] = copy.deepcopy(child_chrsm_1_crs_time)
			next_course_dayslot_population[population_indx] = copy.deepcopy(child_chrsm_1_crs_day)
			next_timeslot_course_population[population_indx] = copy.deepcopy(child_chrsm_1_time_crs)
			next_dayslot_course_population[population_indx] = copy.deepcopy(child_chrsm_1_day_crs)
			population_indx = population_indx + 1

		
		#If infeasible child is produced then
		#rather repairing just discard it. 
		#child 2
			
			
		for crs in 	child_chrsm_2_crs_time.keys() :
			slot = child_chrsm_2_crs_time.get(crs)
			if child_chrsm_2_time_crs.has_key(slot) :
				arr = child_chrsm_2_time_crs.get(slot)
				arr.append(crs)
				child_chrsm_2_time_crs[slot] = arr
			else :
				child_chrsm_2_time_crs[slot] = []
				child_chrsm_2_time_crs[slot] = [crs]
	
		timeslots = child_chrsm_2_time_crs.keys()	
		
		for slot in timeslots:
			crss = child_chrsm_2_time_crs.get(slot)
			#print len(crss) 
			if int(optim_timeslot_or_cost) == 1 and int(len(timeslots)) != max_time_slots :
				select_child_2 = False
				break
			if int(optim_timeslot_or_cost) == 0 and int(len(timeslots)) != min_no_of_timeslots :
				select_child_2 = False
				break
			if len(crss) > 1 :
				for indx1 in range(0,len(crss)) :
					for indx2 in range(0,len(crss)) :
						if int(indx1) != int(indx2) :
							crs1 = int(crss[indx1])
							crs2 = int(crss[indx2])
							#print "2 course conflict "+str(crs1)+"	"+str(crs2)+" val "+str(int(course_conflict_matrix[crs1,crs2]))
							if int(course_conflict_matrix[crs1,crs2]) > 0 :
								select_child_2 = False
								break
					#break
			#break
			if select_child_2 :
				capacity = 0
				for crs in crss :
					capacity = int(capacity) + int(course_student_count_map.get(crs))
		
				if int(capacity) > int(room_capacity):
					select_child_2 = False
					break	
		
		
		#If child 2 produced as a feasible one then populate it to next population
		
		if select_child_2 :
			for timeslot in timeslots :
				slot_key = 0
				if timeslot%every_day_time_slot == 0 :
					slot_key = int(timeslot/every_day_time_slot)
				else :
					slot_key = int(timeslot/every_day_time_slot)+1
			
				if child_chrsm_2_day_crs.has_key(slot_key) :
					child_chrsm_2_day_crs[slot_key].append(child_chrsm_2_time_crs.get(timeslot))
				else :
					child_chrsm_2_day_crs[slot_key] = []
					child_chrsm_2_day_crs[slot_key].append(child_chrsm_2_time_crs.get(timeslot))
			
			next_population_indices[population_indx] = population_indx
			next_course_timeslot_population[population_indx] = copy.deepcopy(child_chrsm_2_crs_time)
			next_course_dayslot_population[population_indx] = copy.deepcopy(child_chrsm_2_crs_day)
			next_timeslot_course_population[population_indx] = copy.deepcopy(child_chrsm_2_time_crs)
			next_dayslot_course_population[population_indx] = copy.deepcopy(child_chrsm_2_day_crs)
			population_indx = population_indx + 1
	
	
	
		
		
		

