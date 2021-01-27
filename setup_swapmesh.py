# I had to split this up into two scripts because of a bug in gmsh.initialize().
# This script creates the setup for swap mesh.
import subprocess

def swap_mesh_setup(l_start=0.02, l_end= 0.12, v_pull=40):  # v_pull in mm/min
    """Create setup for SwapMesh."""
    print('Starting length of crystal:\t', l_start)
    print('Final length of crystal:\t', l_end)
    crys_len = [l_start]
    while crys_len[-1] < l_end:
        crys_len.append(crys_len[-1] + crys_len[-1]/3)
    crys_len[-1] = l_end

    with open('./setup.log', 'w') as f:
        for i in range(len(crys_len)):
            print(f'creating setup for l = {crys_len[i]}')
            subprocess.run(['python', 'setup.py', str(crys_len[i]), f'cz-simple_{i}', str(v_pull)],  stdout=f)
            # model = geometry(crys_len[i], f'cz-simple_{i}')
            # sif(model, v_pull)
            # del model

    with open('./time_vs_length.txt', 'w') as f:
        f.write('time [s]    length [m]\n')
        for length in crys_len:
            time = (length - crys_len[0]) / (v_pull  / 6e4)  # mm / min to m/s
            f.write(f'{time}    {length}\n')


if __name__ == "__main__":
    swap_mesh_setup()
