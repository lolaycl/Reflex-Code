function [T I0_T I1_T I2_T EV_TOTAL HPHI_TOTAL EH_TOTAL] = read_out_file()

filename = ['result.dat'];

imported_data = importdata(filename);

variables_list = imported_data.colheaders;

variables_number = length(variables_list);

number_of_data_per_column = length(imported_data.data);

T = imported_data.data(:,1); % Time vector in micro seconds

EV = imported_data.data(:,2); % Total vertical electric field

EVEL = imported_data.data(:,3); % electrostatic component of the total vertical electric field

EVIN = imported_data.data(:,4); % induction component of the total vertical electric field

EVRD = imported_data.data(:,5); % raditation component of the total vertical electric field

EV_TO = imported_data.data(:,6); % turn-on term of the total vertical electric field

EHO = imported_data.data(:,7); % Total horizontal electric field without Rubinstein correction term

EHF = imported_data.data(:,8); % Total horizontal electric field with Rubinstein correction term

EHEL = imported_data.data(:,9); % electrostatic component of the total horizontal electric field

EHIN = imported_data.data(:,10); % induction component of the total horizontal electric field

EHRAD = imported_data.data(:,11); % raditation component of the total horizontal electric field

EH_TO = imported_data.data(:,12); % turn-on term of the total horizontal electric field

HPHI = imported_data.data(:,13); % Total horizontal magnetic field

HPIND = imported_data.data(:,14); % induction component of the total horizontal magnetic field

HPRAD = imported_data.data(:,15); % raditation component of the total horizontal magnetic field

HPHI_TO = imported_data.data(:,16); % turn-on term of the total horizontal magnetic field

I0_T = imported_data.data(:,19); % Current at the channel base

I1_T = imported_data.data(:,20); % Current at the first point

I2_T = imported_data.data(:,22); % Current at the second point

IFRONT_T_ = imported_data.data(:,23);

EVtower_T_ = imported_data.data(:,28);

EV_TOTAL = EV + EV_TO;

HPHI_TOTAL = HPHI + HPHI_TO;

EH_TOTAL = EHF;