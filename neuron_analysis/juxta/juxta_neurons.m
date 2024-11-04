clear all



cd 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\juxta_BLA\all_principal_neurons\selected'
all_principal_neurons = dir('*selected*');

cd 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\juxta_BLA\LA_neurons\selected'
LA_neurons = dir('*selected*');

cd 'Z:\HajosLab\Dani\Magyar_Daniel\experiments\juxta_BLA\BA_neurons\selected'
BA_neurons = dir('*selected*');

LA_PNs = intersect({all_principal_neurons.name}, {LA_neurons.name})';
BA_PNs = intersect({all_principal_neurons.name}, {BA_neurons.name})';
