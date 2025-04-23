import yaml
import numpy as np

def process_phonon_freq_groupVel(input_data, output_file):
    """
    Extract phonon data1 grouped by band index and save all q-points for each band.
    Args:
        input_data1 (list): List of q-points and their associated data1.
        output_file (str): Path to the output text file.
    """
    # Determine the maximum number of bands from the first q-point
    max_bands = len(input_data[0].get('band', [])) if input_data else 0

    # Open the output file
    with open(output_file, 'w') as out:
        # Write unified header
        out.write("qpoint_index branchID frequency groupVel_x groupVel_y groupVel_z DOS gruneisen_param\n")

        # Write number of band and qpoint
        out.write(f"Max bandID: {max_bands}\n")
        out.write(f"Max qpoint index: {nqpoint}\n")

        # Loop through each band index
        qpoint_index = 1
        for qpoint in input_data:
            # Collect and write data for the current band across all q-points
            for band_index in range(max_bands):
                q_pos = qpoint.get('q-position', [0.0, 0.0, 0.0])
                bands = qpoint.get('band', [])
                
                # Extract data for the current band
                if band_index < len(bands):
                    band = bands[band_index]
                    frequency = band.get('frequency', 0.0)
                    gv = band.get('group_velocity', [0.0, 0.0, 0.0])
                    # Write data for the current band and q-point
                    out.write(
                        f"{qpoint_index} {band_index + 1} "
                        f"{frequency:.10f} {gv[0]:.7f} {gv[1]:.7f} {gv[2]:.7f}\n"
                    )
                    
            qpoint_index = qpoint_index + 1
    print(f"Phonon frequency and group velocity successfully written to {output_file}!")

def process_phonon_dos(input_file,output_file):
    header, data = read_data(output_file, skip_header=3)
    dos_header, dos_data = read_data(input_file, skip_header=1)

    # ensure the number of rows equal
    max_rows = max(data.shape[0], dos_data.shape[0])
    data = pad_data(data, max_rows)
    dos_data = pad_data(dos_data, max_rows)

    dos_col2 = dos_data[:,1]

    # combine header and updated data
    updated_data = np.column_stack((data, dos_col2))

    with open(output_file, 'w') as file:
        file.writelines(header)
        np.savetxt(file,updated_data,fmt='%d %d %.10f %.7f %.7f %.7f %.10f', delimiter=' ') 

    print(f"Phonon DOS sucessfully written to {output_file}!")

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
        # ensure data is 2D before padding
        if data.ndim == 1:
            data = data.reshape(-1,1)

        padding = np.zeros((target_rows-current_rows,data.shape[1]))
        data = np.vstack((data, padding))
    return data

def process_phonon_gruneisen(input_data, output_file):
    """
    Extract phonon data grouped by band index and save all q-points for each band.
    Args:
        input_data (list): List of q-points and their associated data.
        output_file (str): Path to the output text file.
    """
    # Determine the maximum number of bands from the first q-point
    max_bands = len(input_data[0].get('band', [])) if input_data else 0

    # Loop through each band index
    temp=[]
    # Loop through each band index
    for qpoint in input_data:
        # Collect and write data for the current band across all q-points
        for band_index in range(max_bands):
            bands = qpoint.get('band', [])
                
            # Extract data for the current band
            if band_index < len(bands):
                band = bands[band_index]
                gruneisen = band.get('gruneisen', 0.0)
                temp.append(gruneisen)

    # convert temp list to numpy array
    temp = np.array(temp)
    header, data = read_data(output_file, skip_header=3)
    max_rows = max(data.shape[0], temp.shape[0])
    data = pad_data(data, max_rows)
    temp = pad_data(temp, max_rows)

    # combine header and updated data
    updated_data = np.column_stack((data, temp))

    with open(output_file, 'w') as file:
        file.writelines(header)
        np.savetxt(file,updated_data,fmt='%d %d %.10f %.7f %.7f %.7f %.10f %.10f', delimiter=' ')             

    print(f"Phonon Gruneisen parameter successfully written to {output_file}!")

# File paths
freqGroupvel_file = './input/mesh.yaml'
dos_file = './input/total_dos.dat'
gruneisen_file = './input/gruneisen.yaml'
output_file = 'Phonon_Dispersion.dat'

# Main Execution
if __name__ == "__main__":    
    # Load the mesh.yaml file
    with open(freqGroupvel_file, 'r') as f:
        data = yaml.safe_load(f)

    # number of qpoint
    nqpoint = data.get('nqpoint')
    # Extract the phonon data1
    input_data = data.get('phonon', [])
    # Process and save the phonon data grouped by band index
    process_phonon_freq_groupVel(input_data, output_file)

    # load dos file
    process_phonon_dos(dos_file,output_file)

    # Load the YAML file
    with open(gruneisen_file, 'r') as f:
        data = yaml.safe_load(f)

    # Extract the phonon data
    input_data = data.get('phonon', [])
    
    # Process and save the phonon data grouped by band index
    process_phonon_gruneisen(input_data, output_file)
