using WaterLily,ReadVTK
sim = make_sim(...)
# restart the simulation
writer = restart_sim!(sim; fname="file_restart.pvd")
# this acctually append the data to the file used to restart
write!(writer, sim)
# don't forget to close the file
close(writer)