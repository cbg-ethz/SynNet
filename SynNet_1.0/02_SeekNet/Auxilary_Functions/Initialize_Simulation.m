function FoutName = Initialize_Simulation(Constraints)
global Consts
fprintf(['\nNOTE: No "DataPath" found in "' Constraints '", switching to simulaion mode.\n'])
mkdir( Consts.OutFolder);
FoutName = [ Consts.OutFolder '/Simulation at ' strrep(datestr(now), ':', '-') '.txt'];
FoutName = strrep(FoutName, ' ', '_');
disp(['       Simulation Output is Saved at: "' FoutName '"'])
fprintf('\n\n');
% Shuffle MATLAB random stream
s = RandStream('mt19937ar','Seed','shuffle');%s = RandStream('mt19937ar','Seed', 0);
RandStream.setGlobalStream(s);
copyfile(Constraints,[FoutName(1:end-4) '-Used_Constants.log']);
Fout = fopen(FoutName, 'w');
fclose(Fout);
end