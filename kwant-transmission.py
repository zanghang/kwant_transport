import numpy as np
import kwant
import scipy

def read_hami(filename_hami):
    l = []
    atom = []
    atom_single_block = []
    atom_single_block2 = []
    atom_brief = []
    atom_brief2 = []
    atomic_basis = []
    element = []
    hami = []
    atom_number = 0
    block_number = 0
    basis = 0
    #####################################read hamiltonian value###########################################
    with open(filename_hami, 'r') as f:
        for line in f:        
            line = line.split()
            if len(line) >= 4:
                if line[2].isalpha() and len(line[2]) <= 2 : # for identifying the useful line
                    if len(line) == 7 :
                        a, b, c, x1, x2, x3, x4 = line
                        l.append(np.array([float(x1), float(x2), float(x3), float(x4)] ))
                        atom.append([b,c])
                        if int(a) == 1:
                           block_number = block_number + 1
                    elif len(line) == 6 :
                        a, b, c, x1, x2, x3 = line
                        l.append(np.array([float(x1), float(x2), float(x3)] ))
                        atom.append([b,c])
                        if int(a) == 1:
                           block_number = block_number + 1
                    elif len(line) == 5 :
                        a, b, c, x1, x2 = line
                        l.append(np.array([float(x1), float(x2)] ))
                        atom.append([b,c])
                        if int(a) == 1:
                           block_number = block_number + 1
                    else :
                        a, b, c, x1 = line
                        l.append(np.array([float(x1)] ))
                        atom.append([b,c])
                        if int(a) == 1:
                           block_number = block_number + 1
    #print(l)
    
    basis = int(len(atom)/block_number) # get the number of total basis fulction
    
    hami = [[0]*basis for p in range(basis)]
    for i in range(basis):
        for j in range(basis):
            hami[i][j] = float(l[i+int(j/4)*basis][j%4])
    #print(hami)
    #print(l)
    #####################################read hamiltonian value###########################################
    
    ####################get the number of basis function for each atom####################################
    #-------------------------------------------------------
    atom_single_block = [[0] for p in range(basis)]  
    atom_single_block2 = [[0] for p in range(basis)] 
    for i in range(basis):
        atom_single_block[i] = atom[i]                  # get the single block basis , include index and label
        atom_single_block2[i] = atom_single_block[i][1] # get the single block basis , include label , [USEFUL!]
    #print(atom_single_block)  
    #print(atom_single_block2) 
    #-------------------------------------------------------
    for i in atom_single_block2:
        if not i in element:
            element.append(i) # get the unrepeated atom label
    #print(element)
    element_number = len(element) # the number of element
    #-------------------------------------------------------
    for i in atom_single_block: 
        if not i in atom_brief:
            atom_brief.append(i) # get the single block atom , include index and label
    #print(atom_brief) 
    
    atom_number = len(atom_brief) # the number of single block atom
    
    atom_brief2 = [[0] for p in range(atom_number)]
    for i in range(atom_number):
        atom_brief2[i] = atom_brief[i][1] # get the single block atom , include label , [USEFUL!]
    #print(atom_brief2)
    #-------------------------------------------------------    
    dict1 = {}
    for key in atom_brief2:
        dict1[key] = dict1.get(key,0) + 1 # the total atom number for each atom
    #print(dict1)
    dict2 = {}
    for key in atom_single_block2:
        dict2[key] = dict2.get(key,0) + 1 # the total basis number for each atom
    #print(dict2)
    
    key_value = list(dict1.keys())
    key_value2 = list(dict1.values())
    key_value3 = list(dict2.values())
    #print(key_value,key_value2,key_value3)
    
    atomic_basis = [[0]*2 for p in range(element_number)]
    for i in range(element_number):
        atomic_basis[i][0] = key_value[i]
        atomic_basis[i][1] = int(key_value3[i]/key_value2[i])
#    print(atomic_basis)
    ####################get the number of basis function for each atom####################################
    f.close()
    return (hami,atomic_basis,element_number)

def read_coord(filename_hami,filename_coord):
    xyz = []
    atom_list = []
    kwant_family = []
    kwant_family_size = 0
    with open(filename_coord, 'r') as f:
        f.readline()
        f.readline()
        for line in f:
            a, x1, x2, x3 = line.split()
            xyz.append(np.array([float(x1), float(x2), float(x3)]))
            atom_list.append(a)
#    print(atom_list)
#    print(xyz)
    f.close()
    
    hami, atomic_basis, element_number = read_hami(filename_hami) 
#    print(atomic_basis)

############## for repeating the atom xyz to basis xyz i.e. from 88 row xyz to 310 row xyz
    for i in range(len(atom_list)):
        for j in range(element_number):
           if atom_list[i] == atomic_basis[j][0]:
               kwant_family_size += atomic_basis[j][1]
               #kwant_family.append(np.array([float(xyz[i][0]),float(xyz[i][1]),float(xyz[i][2])]))
               #kwant_family.append([float(xyz[i][0]),float(xyz[i][1]),float(xyz[i][2])])
               for k in range(atomic_basis[j][1]):
                   kwant_family.append(tuple([float(xyz[i][0]),float(xyz[i][1]),float(xyz[i][2])]))
    #print (kwant_family_size)
    #kwant_family = [[0.0]*3 for p in range(kwant_family_size)]
    #print (tuple(kwant_family))
#    print (kwant_family)
#    print (hami)
#    for i in range(len(atom_list)):
#        for j in range(element_number):
#           if atom_list[i] == atomic_basis[j][0]:
#               kwant_family[0][]
    f.close()
    return (kwant_family,kwant_family_size,hami)

#===============================================================================================

def kwant_build(device_hami,device_coord,electrode_hami,electrode_coord,threshold,general_lattice,general_device_del,general_device_lead,general_electrode_basis_num):
    basis_sites = []
    
    general_device_del = int(general_device_del)
    general_device_lead = int(general_device_lead)
    general_electrode_basis_num = int(general_electrode_basis_num)
    
    kwant_family, kwant_family_size, hami = read_coord(device_hami, device_coord)
    primitive_vectors = [ (general_lattice,0.0,0.0), (0.0,general_lattice,0.0), (0.0,0.0,general_lattice) ]
    basis_sites = kwant_family
    lat = kwant.lattice.Polyatomic(primitive_vectors, basis_sites)
    syst = kwant.Builder()
    basis = lat.sublattices
    for i in range(general_device_del,kwant_family_size-general_device_del):
        syst[basis[i](0,0,0)] = hami[i][i]# for H1 at 0.0 0.0 0.0

    for i in range(general_device_del,kwant_family_size-general_device_del):
        for j in range(i+1,kwant_family_size-general_device_del) :
            if abs(hami[i][j]) > threshold :
                syst[basis[i](0,0,0),basis[j](0,0,0)] = hami[i][j]
#===================================================================================================
# comment: device top coup with electrode bottom
    lead_family, lead_family_size, lead_hami = read_coord(electrode_hami, electrode_coord)

    sym_up_lead = kwant.TranslationalSymmetry((0.0,0.0,general_lattice))

    up_lead = kwant.Builder(sym_up_lead)

    for i in range(kwant_family_size-general_device_del-general_device_lead,kwant_family_size-general_device_del):
        up_lead[basis[i](0,0,1)] = lead_hami[i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead][i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead]
    
    for i in range(kwant_family_size-general_device_del-general_device_lead,kwant_family_size-general_device_del):
        for j in range(i+1,kwant_family_size-general_device_del):
            if abs(lead_hami[i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead][j+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead]) > threshold :
                up_lead[basis[i](0,0,1),basis[j](0,0,1)] = lead_hami[i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead][j+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead]

    for i in range(kwant_family_size-general_device_del-general_device_lead,kwant_family_size-general_device_del):
        for j in range(kwant_family_size-general_device_del-general_device_lead,kwant_family_size-general_device_del):
            if abs(lead_hami[i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead][j+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)]) > threshold :
                up_lead[basis[i](0,0,1),basis[j](0,0,2)] =  lead_hami[i+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)-general_device_lead][j+general_device_del+general_device_lead-kwant_family_size+int(general_electrode_basis_num/2)] 
    
    syst.attach_lead(up_lead) 
#=======================================================================
# comment: device bottom coup with electrode top
    sym_bottom_lead = kwant.TranslationalSymmetry((0.0,0.0,-general_lattice))

    bottom_lead = kwant.Builder(sym_bottom_lead)

    for i in range(general_device_del,general_device_del+general_device_lead):
        bottom_lead[basis[i](0,0,-1)] = lead_hami[i-general_device_del+int(general_electrode_basis_num/2)][i-general_device_del+int(general_electrode_basis_num/2)]

    for i in range(general_device_del,general_device_del+general_device_lead):
        for j in range(i+1,general_device_del+general_device_lead):
            if abs(lead_hami[i-general_device_del+int(general_electrode_basis_num/2)][j-general_device_del+int(general_electrode_basis_num/2)]) > threshold :
                bottom_lead[basis[i](0,0,-1),basis[j](0,0,-1)] = lead_hami[i-general_device_del+int(general_electrode_basis_num/2)][j-general_device_del+int(general_electrode_basis_num/2)]

    for i in range(general_device_del,general_device_del+general_device_lead):
        for j in range(general_device_del,general_device_del+general_device_lead):
            if abs(lead_hami[i-general_device_del+int(general_electrode_basis_num/2)][j-general_device_del+int(general_electrode_basis_num/2)-general_device_lead]) > threshold :
                bottom_lead[basis[i](0,0,-1),basis[j](0,0,-2)] =  lead_hami[i-general_device_del+int(general_electrode_basis_num/2)][j-general_device_del+int(general_electrode_basis_num/2)-general_device_lead] 

    syst.attach_lead(bottom_lead) 



    def family_colors(site):
#        return 'g' if site.family == basis[0] else 'g'
        return 'b'
    
    kwant.plot(syst,site_size=0.20, site_lw=0.01, hop_lw=0.01,site_color=family_colors,file='zang-test', dpi = 300)

    
    #kwant.plot(syst)
    
    syst = syst.finalized()
    
    #sparse_mat = syst.hamiltonian_submatrix(sparse=False)
    #print(np.matrix(sparse_mat))
    #evs = scipy.sparse.linalg.eigs(sparse_mat,3)[0]
    #print(evs.real)
    energies = []
    data = []
    for ie in range(-80,40):
        energy = ie * 0.01
        smatrix = kwant.smatrix(syst, energy)
        energies.append(energy)
        data.append(smatrix.transmission(1,0))
    outfile = open('transmission.dat','w')
    
    for i in range(len(energies)):
        outfile.write('{0:40.30f} {1:40.30f} \n' .format(energies[i],data[i]))
    outfile.close()


kwant_build(
    device_hami='../cp2k-device-2Au/TEST-zang-device-1_0_502.Log', 
    device_coord='../cp2k-device-2Au/device.xyz',
    electrode_hami='../cp2k-electrode/TEST-zang-electrode-1_0_502.Log',
    electrode_coord='../cp2k-electrode/electrode.xyz',
    threshold = 0.0,
    general_lattice = 12.0,
    general_device_del = 8*9,
    general_device_lead = 4*9,
    general_electrode_basis_num = 24*9,
    )
