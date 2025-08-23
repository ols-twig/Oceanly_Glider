basic minimum file structure for processing slocum glider files
- cac						-- where all the cache files for the glider should be stored
- raw_dbd_ebd				-- where to put all the *.dbd and *.ebd files downloaded from the glider 
- rawnc						-- where all the processed files will be stored after running 'run_program_xxx.py'
- deploymentRealtime.yml	-- change this to contain all the mets information for the current project
- unit_1104_sensors.txt		-- a list of all the masterdata parameters to be processed into the *.nc's