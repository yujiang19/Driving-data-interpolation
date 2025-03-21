# Driving-data-interpolation


A solver for interpolation of low-frequency driving data. The interpolation algorithm based on driver behavior characteristics is explained in detail in the paper:

[1] Electric vehicle driving data reconstruction and validation: A method combining behavioral characteristics and multi-perspective insights

# **Overview**

Big data platforms provide strong support and network advantages for innovating and optimizing energy-saving technologies for electric vehicles (EVs). However, it faces the dual challenges of low-resolution data and processing complexity. We proposed a comprehensive method for data reconstruction and validation incorporating behavioral characteristics and multi-perspective insights to address these challenges. The proposed method introduces driving characteristics extraction, pattern recognition, and reconstruction module design to achieve data interpolation tasks.

# Main file description

This respository contains the following files:

**example_inpdata.mat -**  A Matlab data file that contains the inputs for the solver

**main.m -**  Main function

**augmentf.m -**  Driving data interpolation solver as Matlab function file, including vehicle speed interpolation, current interpolation, and corresponding baseline model

       valdationpolt.m is a class function that includes the following methods
       function obj=augmentf(D)
            obj.D=D; % 
            obj.initial; % initialization
            obj.vagment; % Our imputation model
            obj.linear; % Linear interpolation model
            obj.spline; % Spline interpolation model
            obj.speed_prcess; % BBI and markov-based interpolation model
            obj.SOCprocess; % Battery data interpolation model
        end

**valdationpolt.m –**  A functions for model validation and result visualization

       valdationpolt.m is a class function that includes the following methods：

**obj.calfea -**  Calculate characteristic parameter 

**obj.vplot -**   Vehicle speed calibration visualization

**obj.socplot -**   Battery data visualization

**obj.vaplot -**   Calculate velocity-acceleration distribution

**obj.energyplot -**   Calculate energy consumption

Other functions should also include the driving speed prediction function and energy consumption calculation function for interpolation model verification. Unfortunately, since the dataset used for training the above functions is proprietary, I cannot disclose the code of the entire model. I will disclose this part of the code after subsequent modifications.

# Input data

Due to confidentiality requirements for the original driving data, we cannot make the complete dataset public. Therefore, we provide sample data that contains the key variables required for the model to run, including vehicle speed, acceleration, driving distance, accelerator and brake pedal opening, and battery data. We hope that this example can help users understand the basic structure and format requirements of the data required in this research model, and provide a reference for users to verify it with their own data.

The data_trip1s.mat file in the example_data folder is data at a sampling frequency of 1HZ, and the data_trip10s.mat file is data at a sampling frequency of 0.1HZ.

       data_trip10s.mat包含了以下数据：

**Column6 -**   A Nx1 array of vehicle driving speed in m/s

**Column8 -**   A Nx1 array of the battery voltage in V

**Column9 -**   A Nx1 array of battery current in A

**Column10 -**   A Nx1 array of SOC

**Column13 -**   A Nx1 array of pedal opening

**Column16 -**   A Nx1 array of vehicle acceleration in m/s2

**Column17 -**   A Nx1 array of driving distance in m

The Current fitting function.mat file in the similar driver data folder is the fitting function based on high-frequency driver data and the driver_limit.mat file is the similar driver constraint information.

# Results

1、Speed reconstruction results under different methods

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/1.png" width="1500px" alt="fig1">
</div>

2、Comparison of the error between the reconstructed speed and the original speed in characteristic parameters using different methods. (The horizontal axis labels from left to right are: 1(maximum speed), 2(average speed), 3(average driving speed), 4(standard deviation of speed), 5(average acceleration), 6(standard deviation of acceleration), 7(average positive acceleration), 8(maximum positive acceleration), 9(maximum deceleration), 10(average deceleration), 11(driving distance), 12(percentage of driving time under positive acceleration), 13(percentage of driving time under deceleration), 14(percentage of idle time).)

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/2.png" width="500px" alt="fig2">
</div>

3、Resampling results of current and voltage

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/3.png" width="800px" alt="fig3">
</div>

4、The kernel density function curves of acceleration across the entire speed range.

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/4.png" width="500px" alt="fig4">
</div>

5、SOC resampling results.

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/5.png" width="800px" alt="fig5">
</div>

6、Comparison between simulated energy consumption and measurements.

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/6.png" width="500px" alt="fig6">
</div>

7、Error box plot of battery calibration model results.

<div align=center>
<img src="https://github.com/yujiang19/pic/blob/main/7.png" width="500px" alt="fig7">
</div>

# Disclaimer

This code is provided as freeware and is intended solely for non-commercial, educational, and research purposes. It must not be used for any commercial purposes without prior authorization from the code developer. Any use for commercial purposes without such authorization will render you and the users responsible for any resultant liabilities, and the code developer and the platform will not be held responsible for any consequences arising therefrom. Users assume all risks associated with the use of this code. The developer and associated platforms disclaim any liability for special, incidental, direct, or indirect damages arising out of or in connection with the use or inability to use the code. This includes, but is not limited to, any loss of data or property, and any resulting or related liabilities to the user or any third parties. By downloading or using this code, you signify your agreement to these terms.

# Contact Us

Please contact us if you need further technical support. If you have any trouble with this repo, feel free to contact us by e-mail. We'll try to resolve the issue as soon as possible! Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Email contact: [Yu Jiang](15615630556@163.com)

Ph.D. Candidate, Jilin University.
