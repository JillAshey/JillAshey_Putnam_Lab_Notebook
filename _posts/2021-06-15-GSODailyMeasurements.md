---
layout: post
title: Daily measurements 
date: '2021-06-15'
categories: Protocols
tags: [Protocols]
projects: GSO Astrangia 
---

## GSO *Astrangia* Experiment: Daily Measurements 

**About**: This post outlines daily tasks to be completed at GSO for J. Ashey's *Astrangia* experiment. It is adapted from [protocols](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/GSO_Wetlab_Protocols/GSO_Wetlab_Protocols.md) created by DMB and SJG. More information about OA experiments and carbonate chemistry is also in these protocols. 

### Equipment needed 

| Instrument  | Measurement | Measurement frequency
| ------------- | ------------- | ------------- | 
| Traceable Digital Thermometer | Temperature (°C) | Daily
| Orion Star A325 Thermoscientific pH/conductivity meter | pH (mV) and conductivity (psu) | Daily
| Apogee underwater quantum meter (MQ-510) | Light (µmol m-2s-1) | 2x week
| HOBO pendant MX2202 data logger | Temperature (°C) and light (lux) | Read out weekly

### Before getting started...

Think critically about your experiment

- What is the duration of your experiment?
- How often are you sampling/measuring responses?
- What is the nature of your treatment(s) and experimental design?

Record all data in lab notebook or Astrangia binder. Be sure to write **initials and the date for each entry**.

## Daily tasks

### I. Upon arrival 

1.**Check the tanks and equipment visually**. Make sure that...

- All lights are on; if one is turned off, unplug it, wait one minute and plug it back in. 
- All chiller tubing is properly placed in ambient tanks (i.e., not falling out) and one of the pumps in the ambient tank is connected to the tubing. 
- All equipment is plugged in, either to the outlets on the back of the tanks, the outlets against the back wall, or the Apex system. 

2.**Check the Apex Fusion app** on your phone or computer. Make sure these data match the desired conditions in each tank.


### II. Calibrations 

#### Calibration overview 

| Calibration  | Instruments | Measurement frequency | Solutions  
| ------------- | ------------- | ------------- | ------------- | 
pH | Orion Star A325 Thermoscientific pH/conductivity meter (pH probe); TEMP PROBE | 2x week | Dickson Lab Tris standard
Conductivity | Orion Star A325 Thermoscientific pH/conductivity meter (cond. probe) | 2x week | Conductivity standard, 50,000 µS/cm 
Apex temperature probes | Temp PROBE | 1x week | NA

#### pH Tris calibrations

**About**: We need to measure tris standard in order to convert our pH (mV) discrete measurements to total scale via Dickson's [Guide to Best Practices](https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/Handbook_2007.html). The goal of the Tris calibration is to measure pH (mV) across a range of temperatures ± 5-10 Celcius than currently experienced in the tanks. The temperature values and corresponding pH values must have a linear relationship and an R^2 value equal to or above 0.98 in order to calculate the pH values for the daily measurements from the calibration curve. Do this calibration 2x a week.

1. Fill a 50 mL falcon tube with ~30-35 mL of certified Tris standard solution from the Dickson lab. 
	- Replace the Tris solution 2x a week.
2. Fill a 1L tripour beaker with water and put a frozen water bottle in it. Position the falcon tube with the Tris solution in the tripour beaker next to the frozen bottle to cool the Tris solution down to the desired temperature. If there are not frozen bottles or ice available, put the falcon tube in a refridgerator or freezer. 
	-  Take note of your desired temperature. This value will change based on the range of temperatures you expect to be measuring in your tanks. You want the temperature of the Tris solution to decreased to ~5°C below the lowest temperature in your experiment at that point in time. 
3. With the temperature probe, measure the Tris solution temperature. Swirl the temperature probe in the solution, as it measures to ensure the solution is well mixed.
4. Once the solution reaches the desired temperature, turn on the Orion pH meter and put the pH probe in the falcon tube with the temperature probe. Swirl both probes in the solution as it measures to ensure the solution is well mixed. 
5. Click 'Measure' on the Orion meeter when you are ready to take your first measurement. 
6. Wait for mV value to settle (will say 'ready'), and record the pH and temperature in the lab notebook.
7. Warm the tube of Tris solution by gently holding the tube in your hand. Be careful to not hold the container for too long or too aggressively, as the temperature will increase too fast for calibration measurements. 
8. Record temperature and pH measurements across the desired range, with increments of about ~0.5-1°C. Read the temperature value as soon as the pH meter states 'ready' instead of 'stabilizing'.
9. Start a new CSV file (using the header below) named today's date (yyyymmdd.csv) and input the pH and temperature calibration measurements. Save the data to your pH_Tris data folder. 
	- Use these headers in the CSV file:

| mVTris | TTemp |
| ------------- | ------------- |


#### Conductivity calibration

#### Apex temperature probe calibration

1. Turn the temperature probe on and put the metal tip in the tank. 
2. Wait for the probe to settle within 0.1°C and record the temperature in the lab notebook (can do this during daily measurements).
3. Open the Apex Fusion app on your phone or computer. Navigate to the dashboard displaying temperatures for all tanks. 
4. Click the wheel (for Settings) in the top right corner of the box displaying temperature probe name and current tank temperature for desired tank. It will take you to the 'Temperature Probe Configuration' page.
5. Click on the orange 'Calibration' button at the bottom of the page.
6. A smaller window will pop up. When you are ready to start the calibration, click 'Calibrate'.
7. Enter the temperature obtained from the temperature probe. Click 'Next'. 
8. The window should now read 'Calibration Complete'. Click 'Done'.
9. In the top right cornor of the screen, click the button with a cloud and an arrow within it. This button will update the Apex.
10. Repeat steps 3 - 9 for all replicate tanks and head tanks. 
11. After calibrations for all tanks are completed, turn off the temperature meter, rinse the probe with DI water, and dry with a kimwipe. 


### III. Discrete measurements 

Use the Astrangia binder with pre-printed data sheets to record daily measurements. The data sheet headers should look like: 

| Date | Time | Tank | Temperature | pH (mV) | Salinity | Light | Initials |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | 

#### Temperature measurements

1. Wipe off the metal portion of the probe with DI water using a kimwipe
2. Press the "ON/OFF" button to turn the probe on.
3. Place the tip/metal portion of the probe in the tank.
4. When the probe settles within 0.1°C, record the temperature in the lab notebook
5. Repeat steps 3 & 4 for all replicate tanks and head tanks. NOTE: record all temperature measurements in Celsius. Press the F/C button to reset as necessary
6. Press the "ON/OFF" button to turn the probe off (and complete the daily temperature probe care checklist below).
7. Rinse metal probe with DI water and dry with a kimwipe. 
8. Wipe off meter box with a kimwip with 70% ethanol to remove water or salt. 

#### pH measurements 
1. Press the power on the Orionn meter to turn the probe on
2. Rinse the tip of the probe with DI water and dry with Kimwipe. BE VERY CAREFUL, this is a glass electrode probe
3. Toggle the "channel" button (F3) to make sure you are viewing pH in mV
4. Submerge the tip of the probe in each tank to take recording. Do not submerge the probes too far to avoid corrosion or short-circuiting connections.
5. When the probe settles (will say "ready") record the pH in the lab notebook
6. Repeat steps 4 & 5 for each replicate tank and head tank. NOTE: record all temperature measurements in mV. Press the "mode" button to reset as necessary
7. Press the power button to turn the probe off (and complete the daily salinity probe care checklist below).
8. Rinse pH probe with DI water and dry with a kimwipe. Be very careful when drying the glass electrode!
9. Store pH probe in the pH electrode storage solution with the electrode completely submerged in the solution (no bubbles around the electrode).
10. Wipe off Orion meter with a kimwip with 70% ethanol to remove water or salt. 

#### Conductivity (salinity) measurements 

1. Press the power button to turn the probe on
2. Toggle the "channel" button (F3) to make sure you are viewing salinity (conductivity) in psu
3. Place the probe in the tank. Make sure the whole conductivity cell is submerged in the water The conductivity probe needs to be positioned below water just enough that the cell/opening is submerged.
4. When the probe settles (will say "ready") record the salinity in the lab notebook
5. Repeat steps 3 & 4 for all replicate tanks and head tanks. NOTE: record all salinity measurements in psu. Press the "mode" button to reset as necessary
6. Press the power button to turn the probe off.
7. Rinse salinity probe with DI water and dry with a kimwipe.
8. Store the salinity probe in DI water with the conductivity cell covered in DI water
9. Wipe off Orion meter with a kimwip with 70% ethanol to remove water or salt.

#### Light measurements 

1. Press the power button to turn on the Apogee meter and remove the cap on the sensor.
2. Using the extension stick, put the sensor in the middle of the coral rack. 
3. The sensor is very sensitive and will fluctuate rapidly. Take the average of the values you are seeing and record the light values in the lab notebook. 
4. Repeat steps 3 & 4 for all replicate tanks and head tanks.
5. Press the power button to turn the Apogee meter off. 
6. Rinse sensor with DI water and dry with a kimwipe. Put cap back on sensor. 
7.  Wipe off Apogee meter with a kimwip with 70% ethanol to remove water or salt.


### IV. Logger readout

1. Remove logger from tank and rinse it off in the sink. Scrub the algae off if needed.
2. Open the HoboConnect app on your phone. You should be on the Devices screen and there should be a radiating bulls-eye in the middle of the screen.
3. To connect the logger with the app, press hard with your thumb on the middle of the pendant. If the logger and app connected, a red light and a green light will flash on the pendant. The logger icon and ID will also show up on the app. 
4. Click on the logger icon. 
5. This will take you to the configuration details of that particular logger. At the bottom of the screen, click the icon that has a square with an arrow pointing out to read the data from the logger. The screen will say 'Readout Complete' when this is finished. 
6. Navigate back to the Devices screen. At the bottom of the screen, next to 'Devices', click on the HOBO files icon. This will take you to your data files 
7. At the top right hand of the screen, click the folder icon. This will give you a list of the data files most recently downloaded. Click on the file that corresponds to the serial number of the specific logger. 
8. Youl'll now see temperature and light data from the logger on a plot. To export the data, click the three dots in the upper right hand corner of the screen. Then click the icon with the square and the arrow pointing out of the square. 
9. Click 'Export to CSV'. 
10. Once the exporting is completed, click 'Share'. You can either airdrop it to your computer or email it to yourself. Once the file is on your computer, rename it with the tank number and date. 


### V. Data logging and backup

1. Take pictures of lab notebook and Astrangia binder and upload them to the Google Drive. Make sure date, initials, and page number is visible in every picture. 
2. Update [Daily CSV file](https://github.com/JillAshey/Astrangia_repo/blob/main/data/DailyMeasurements.csv) with data you just collected and save the file. 
3. With the updated CSV, run the Daily R Markdown [script](https://github.com/JillAshey/Astrangia_repo/blob/main/scripts/Daily.Rmd). This will generate new daily measurement plots that can be uploaded to Github. 


### VI. General check & reminders

1. Check flow in all tanks. If flow is too low or high, adjust using the orange knob. Keep in mind that adjusting the flow in one tank will (slightly) affect the flow going to the other tanks. Depending on the incoming water temperature, the flow rate will affect the tank temperatures. If the incoming water temperature is affecting the tank temperatures, adjust flow rates accordingly.
2. Using the suction cups, put the heaters back on the side of the tank if they have fallen off. 
3. Wipe off the light loggers with thumb so they can accurately measure light. 
4. Before leaving for the day, check the Apex Fusion app to ensure temperatures are still in your desired range. 