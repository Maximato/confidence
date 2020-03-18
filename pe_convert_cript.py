import os


data_dir = "data_concensuses"
files = os.listdir(data_dir)
for file in files:
    command1 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}1_mut90.pe -ml c90 -cf 50 -ct 2000 -fmt pe"
    command2 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}2_mut90.pe -ml c90 -cf 1900 -ct 3850 -fmt pe"
    command3 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}3_mut90.pe -ml c90 -cf 3800 -ct 5750 -fmt pe"
    command4 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}4_mut90.pe -ml c90 -cf 5700 -ct 7650 -fmt pe"
    command5 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}5_mut90.pe -ml c90 -cf 7600 -ct 9550 -fmt pe"
    command6 = f"python mutate.py -i {os.path.join(data_dir, file)} -o {file.split('.')[0]}6_mut90.pe -ml c90 -cf 9500 -ct 11400 -fmt pe"
    os.system(command1)
    os.system(command2)
    os.system(command3)
    os.system(command4)
    os.system(command5)
    os.system(command6)
