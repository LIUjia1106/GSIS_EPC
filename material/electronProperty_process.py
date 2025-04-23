import numpy as np

def process_electron_energy_groupvelocity(input_file, output_file):
    # load data
    data = np.loadtxt(input_file)

    # extract columns
    kpoint_index = data[:, 0].astype(int)
    bandID = data[:, 1].astype(int)
    other_columns = data[:, 2:]

    # find maximum bandID and kpoint_index
    max_bandID = np.max(bandID)
    max_kpoint_index = np.max(kpoint_index)

    # sort with bandID
    # sorted_indices = np.lexsort((kpoint_index, bandID))
    # sorted_data = data[sorted_indices]

    # adjust column order, 1st column is bandID, 2nd column is kpoint_index
    # sorted_data[:, [0, 1]] = sorted_data[:, [1, 0]]
    sorted_data = data

    # write data into file
    with open(output_file, 'w') as f:
        # file head
        f.write("kpoint_index bandID energy groupVel_x groupVel_y groupVel_z DOS fermiEnergy\n")
        f.write(f"Max kpoint index: {max_kpoint_index}\n")
        f.write(f"Max bandID: {max_bandID}\n")
    
        # save data
        np.savetxt(f, sorted_data, fmt='%d %d %.6f %.6f %.6f %.6f', delimiter=' ')

    print(f"Electron energy and group velocity sucessfully written to {output_file}!")

def process_electron_dos_fermienergy(input_file, output_file):
    header, data = read_data(output_file, skip_header=3)
    dos_header, dos_data = read_data(input_file, skip_header=0)

    # ensure the number of rows equal
    max_rows = max(data.shape[0], dos_data.shape[0])
    data = pad_data(data, max_rows)
    dos_data = pad_data(dos_data, max_rows)

    dos_col2 = dos_data[:,1]
    dos_col3 = dos_data[:,2]
    updated_data = np.column_stack((data, dos_col2, dos_col3))


    with open(output_file, 'w') as file:
        file.writelines(header)
        np.savetxt(file,updated_data,fmt='%d %d %.6f %.6f %.6f %.6f %.15f %.6f', delimiter=' ') 

    print(f"Electron DOS and Fermi energy sucessfully written to {output_file}!")

# Extract file data
def read_data(file_path, skip_header):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        header = lines[:skip_header]
        data = np.loadtxt(lines[skip_header:])
    return header, data

# Fill blank data
def pad_data(data, target_rows):
    current_rows = data.shape[0]
    if current_rows < target_rows:
        print("# of rows not match! Autofilled with 0!")
        padding = np.zeros((target_rows-current_rows,data.shape[1]))
        data = np.vstack((data, padding))
    return data

# file path
energyGroupvel_file = './input/electron_gv_200.000007464979.txt'
dosFermi_file = './input/electron_gv_200.000007464979.txt'
output_file = 'Electron_Band.dat'

# Main Execution
if __name__ == "__main__":
    # process and save the electron energy and group velocity sorted by band index
    process_electron_energy_groupvelocity(energyGroupvel_file, output_file)

    # add electron dos and Fermi energy
    process_electron_dos_fermienergy(dosFermi_file,output_file)