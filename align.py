import sys
import getopt
import numpy as np
import rmsd
import shutil
import os


def divide_name(filename):
    #home-made os.path.splitext, for it can't handle "name.a.b.c" properly
    basename = os.path.basename(filename)
    parts = basename.split(".") #split return >= 1 length list
    if len(parts) == 1:
        return parts[0], ""
    else:
        return parts[0], "."+".".join(parts[1:]) 
def get_cell_name(filename):
    cell_name, _ = divide_name(filename)
    return cell_name

THRESHOLD=2
def do_pick(pairs:list,dv:list):
    #pick out structure dv bigger than THRESHOLD, pairs: all combinations in order, dv: median_devariations
    ##get all structure name(represent by int index) 
    problematic = []
    [problematic.extend(i) for i in pairs]
    problematic = set(problematic)
    ##do the pick
    for i,v in enumerate(dv):
        if v < THRESHOLD:
            problematic -= pairs[i]
    #print(problematic)
    ##get index of good pairs
    good_pairs = [i for i,j in enumerate(pairs) if len(j.intersection(problematic)) == 0]
    return problematic, good_pairs
def align_main(args):
    filenames, output_dir, good_dir, bad_dir = args.filenames, args.output_dir, args.good_dir, args.bad_dir
    os.makedirs(good_dir,exist_ok=True)
    os.makedirs(bad_dir, exist_ok=True)    
    # ------------------------load 3dg files--------------------------
    input_data = []
    num_structures = len(filenames)
    if num_structures < 2:
        sys.stderr.write("[E::" + __name__ + "] at least 2 structures are required\n")
        return 1
    counter = 0
    for input_filename in filenames:
        sys.stderr.write("[M::" + __name__ + "] reading 3dg file " + str(counter) + ": " + input_filename + "\n")
        input_data.append({})
        for input_file_line in open(input_filename, "rb"):
            input_file_line_data = input_file_line.strip().split()
            #store in the newly added empty dictionary
            input_data[-1][(input_file_line_data[0], int(input_file_line_data[1]))] = [float(input_file_line_data[2]),float(input_file_line_data[3]),float(input_file_line_data[4])]
        counter += 1
    # --------------------------------find common particles--------------------------------
    # find common particles
    common_loci = set(input_data[0])
    ##get all common_loci
    for input_structure in input_data[1:]:
        common_loci = common_loci.intersection(set(input_structure))
    num_loci = len(common_loci)
    common_loci = list(common_loci)
    ##select data subset for each structure according to common_loci
    common_data = []
    for input_structure in input_data:
        common_data.append([])
        for common_locus in common_loci:
            common_data[-1].append(input_structure[common_locus])
    
    sys.stderr.write("[M::" + __name__ + "] found " + str(num_loci) + " common particles\n")
    
    #common_data = np.array(common_data)
    # --------------------subtract centroid---------------------
    common_data = np.array(common_data)
    centroid_data = []
    #normalize to centroid for each structure
    for i in range(num_structures):
        common_data[i] = np.array(common_data[i])
        centroid_pos = rmsd.centroid(common_data[i])
        common_data[i] -= centroid_pos
        centroid_data.append(centroid_pos)
    sys.stderr.write("[M::" + __name__ + "] found centroids for " + str(num_structures) + " structures\n")
    # ------------calculate pairwise deviation and rotate------------
    #for deviations: locus:sturcture
    deviations = np.empty((num_loci, 0), float)
    #for pick up
    median_deviations = []
    dv_pairs = []
    for i in range(num_structures):
        for j in range(num_structures):
            if j == i:
                continue
            # mirror image if needed
            mirror_factor = 1.0
            if rmsd.kabsch_rmsd(common_data[i], common_data[j]) > rmsd.kabsch_rmsd(common_data[i], -1.0 * common_data[j]):
                mirror_factor = -1.0
            # calculate deviation
            rotation_matrix = rmsd.kabsch(mirror_factor * common_data[j], common_data[i])
            if j > i:
                #print("matrix",np.dot(mirror_factor * common_data[j], rotation_matrix) - common_data[i])
                deviation = np.linalg.norm(np.dot(mirror_factor * common_data[j], rotation_matrix) - common_data[i], axis = 1).T
                #print("deviation",deviation)
                deviations = np.c_[deviations, deviation]
                #print("deviations",deviations)
                dv_pairs.append( set([i,j]) )
                median_deviations.append(np.median(deviation))
                print("[M::" + __name__ + "] median deviation between file " + str(i) + " and file " + str(j) + ": " + str(np.median(deviation)) + "\n")
            
            # rotate
            if output_dir is not None:
                # rotate j to align with i
                sys.stderr.write("[M::" + __name__ + "] aligning file " + str(j) + " to file " + str(i) + "\n")
                cell_name = get_cell_name(filenames[0])
                aligned_filename = os.path.join(output_dir, cell_name + "." + str(j) + "_to_" + str(i) + ".3dg")
                aligned_file = open(aligned_filename, "w")
                for input_locus in input_data[j]:
                    aligned_pos = np.dot((np.array(input_data[j][input_locus]) - centroid_data[j]) * mirror_factor, rotation_matrix) + centroid_data[i]
                    aligned_file.write("\t".join( [str(input_locus[0]), str(input_locus[1]), str(aligned_pos[0]), str(aligned_pos[1]), str(aligned_pos[2])] ) + "\n")
                aligned_file.close()
    

    # ------------exclude structure deviation bigger than threthold------------
    problematic, good_pairs = do_pick(dv_pairs, median_deviations)
    if len(problematic) > 0:
        exclude_files = [filenames[i] for i in problematic]
        good_files = [name for name in filenames if name not in exclude_files ]
        for name in good_files:
            shutil.copy(name, os.path.join(good_dir, os.path.basename(name)) )
        for name in exclude_files:
            shutil.copy(name, os.path.join(bad_dir, os.path.basename(name)) )
        sys.stderr.write("M:: %s: copy good files %s\n" % (__name__, str(good_files)) )
        if len(good_files) >= 3:
            print("Exclude: ", str(exclude_files), "from rmsd calculation.")
            print("Using: ", str(good_files), ".") 
            good_deviations = [deviations[:,i] for i in good_pairs]
            good_deviations = np.array(good_deviations).T
            #print("deviations", deviations.shape)
            #print("good_deviations", good_deviations.shape)
            deviations = good_deviations
        else:
            print("\n\n\nBroken structure( RMSD high for MOST pairs. Unable to pick 3 from within.).\n\n\n")
    else:
        print("All structure good.")
        for name in filenames:
            shutil.copy(name, os.path.join(good_dir, os.path.basename(name)) )   
    #print(dv_pairs)
    #print(median_deviations)

    # ------------calculate rmsds------------
    ## rms between cells
    ## then rms and median between locus
    RMS = lambda dv, i=0 : np.sqrt( (dv**2).mean(axis=i) )
    Median = lambda dv, i=0 : np.median( dv, axis=i )
    median_rmsd = Median( RMS(deviations,1) )
    rms_rmsd = RMS( RMS(deviations,1) )
    
    # ------------print and log------------ 
    print("[M::" + __name__ + "] RMS RMSD: " + str(rms_rmsd))
    print("[M::" + __name__ + "] median RMSD: " + str(median_rmsd))
    sys.stderr.write("[M::" + __name__ + "] writing output\n")
    
    '''
    #print rmsd per loci
    for i in range(num_loci):
        sys.stdout.write("\t".join(map(str, [common_loci[i][0], common_loci[i][1], rmsds[i]])) + "\n")
    '''
     

    