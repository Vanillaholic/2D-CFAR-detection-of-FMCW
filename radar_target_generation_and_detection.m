clear all;
clc;

%% 参数设置
max_range = 200;           %最大距离
range_resolution = 1;      %距离分辨率
max_velocity = 70;         %最大速度
c = 1500;                  %声速


%% 目标的初始位置和速度
target_range = 110;                %目标初始距离（最大200m）
target_velocity = -10;             %-10 m/s（负号表示目标正向雷达靠近）-70~70
 
%% FMCW 波形生成

B = c /(2*range_resolution);       %带宽
Tchirp= (5.5*2*max_range)/c;       %计算 Chirp 时间（5.5 倍往返时间）
slope = B/Tchirp;                  %FMCW斜率
disp(slope);                       % 显示斜率


fc= 1.5e3;                         %载波频率
                                                          
% 一组 Chirp 信号的数量（建议使用 2 的幂次方以便 FFT 计算多普勒频谱）
Nd=128;                            % 多普勒维度（即 Chirp 数量）

%每个chirp的采样数
Nr=1024;                           % 距离维度（即每个 Chirp 的采样点数）

% 时间向量，覆盖所有 Chirp 和每个 Chirp 内的样本
t=linspace(0,Nd*Tchirp,Nr*Nd);     % 生成总时间轴

%% 初始化信号向量
Tx=zeros(1,length(t));            % 发射信号
Rx=zeros(1,length(t));            % 接收信号
Mix = zeros(1,length(t));         % 混频信号（打拍频信号）

% 初始化目标位置和时间延迟的向量
r_t=zeros(1,length(t));           % 记录每个时间点的目标位置
td=zeros(1,length(t));            % 记录每个时间点的时间延迟


%% 信号生成与移动目标仿真
% 在时间轴上运行的场景

for i=1:length(t)         
    
    % *%TODO* :
    % 根据时间步长更新目标的距离 (匀速运动)
    % r(t) = r0 + v*t
    r_t(i) = target_range+ target_velocity*t(i); % 目标当前距离
    td(i) = 2*r_t(i)/c;                          % 计算信号往返时间延迟
    % TODO:
    % 为每个时间点生成发射和接收信号
    % 发射信号 Tx：FMCW 信号的基本公式
    Tx(i) = cos(2*pi*(fc*t(i) + (slope*t(i)^2)/2));

    % 接收信号 Rx：考虑了往返延迟的 FMCW 信号
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i)) + (slope*(t(i)-td(i))^2)/2));
    
    % TODO:
    % 生成混频信号（打拍频信号）
    % 将发射信号和接收信号逐元素相乘
    Mix(i) = Tx(i).*Rx(i);
    
end

%% 距离测量 (Range Measurement)


% 将混频信号（Mix）重塑为 Nr x Nd 的矩阵
% Nr：距离维度（每个 Chirp 的采样点数）
% Nd：多普勒维度（Chirp 数量）
sig = reshape(Mix, [Nr, Nd]);


% 在距离维度 (Nr) 上对打拍频信号进行 FFT 变换
sig_fft1 = fft(sig, Nr); 

% 对 FFT 结果进行归一化处理
sig_fft1 = sig_fft1 ./ Nr;


% 取 FFT 结果的幅值（只关注频谱强度）
sig_fft1 = abs(sig_fft1);


% FFT 结果是双边谱，只保留单边谱（正频率部分）
sig_fft1 = sig_fft1(1:(Nr/2));

% 绘制距离维度的 FFT 结果
figure('Name','Range from First FFT')

% 绘制 1D FFT 输出，展示距离信息
plot(sig_fft1, "LineWidth",2);
grid on;axis ([0 200 0 0.5]);
xlabel('range');ylabel('FFT output');title('1D FFT');



%% 距离-多普勒响应 (Range-Doppler Response)
% 该部分实现了 2D FFT，用于生成距离-多普勒图 (RDM)。
% 后续可以基于该 RDM 进行 CFAR 检测。


%二维 FFT 的输出是一个包含距离和多普勒频移响应的图像，其坐标轴单位最初是FFT bin（频率单元格）。
% 因此，为了直观反映目标的实际距离和速度，
% 需要根据最大取值将坐标轴从 bin 索引转换为实际的距离 (m) 和多普勒速度 (m/s)。


% 将混频信号重塑为 Nr × Nd 矩阵
Mix=reshape(Mix,[Nr,Nd]);

% 对打拍频信号执行二维 FFT（分别在距离和多普勒维度上）
sig_fft2 = fft2(Mix,Nr,Nd);


% 仅保留距离维度的单边谱
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);

sig_fft2 = fftshift(sig_fft2);% 对多普勒维度做 fftshift，使零频分量居中
RDM = abs(sig_fft2);        % 计算幅度并转换为 dB 值
RDM = 10*log10(RDM) ;

% ---------------------
% 绘制距离-多普勒图 (RDM)
% ---------------------

% 定义多普勒和距离坐标轴
doppler_axis = linspace(-100,100,Nd);             % 多普勒频移轴（-100 ~ 100 m/s）
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);% 距离轴

figure;
surf(doppler_axis,range_axis,RDM);
xlabel('doppler');ylabel('range');zlabel('RDM');title('2D FFT');


%% CFAR implementation

%在RDM上窗口滑动

% 选择训练单元 (Training Cells) 的数量
Tr = 8; % 距离维度上的训练单元数量
Td = 4; % 多普勒维度上的训练单元数量

% 选择保护单元 (Guard Cells) 的数量，围绕 CUT 以避免信号泄漏
Gr = 4; % 距离维度上的保护单元数量
Gd = 2; % 多普勒维度上的保护单元数量

% *%TODO* :
% 设置SNR偏移量 (以dB为单位)，用于调整检测灵敏度
offset = 1.4;
% *%TODO* :
%初始化一个向量来存储训练单元每个 iteration 的噪声级
noise_level = zeros(1,1);


% *%TODO* :
% 设计一个循环，使 “CUT” 在距离 - 多普勒图上滑动，同时在边缘处为训练单元和保护单元留出边界。
% 对于每次迭代，将所有训练单元内的信号电平进行求和。为了求和，需使用 db2pow 函数将数值从对数形式转换为线性形式。
% 对所使用的所有训练单元的求和值求平均。求平均后，再使用 pow2db 函数将其转换回对数形式。
% 进一步给它加上偏移量以确定阈值。接下来，将 “CUT” 下的信号与该阈值进行比较。
% 如果 “CUT” 的电平 > 阈值，则赋予其值为 1，否则将其设为 0。


   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR

RDM = RDM/max(max(RDM));
 %滑动窗口遍历整个RDM
for i = Tr+Gr+1 : Nr/2-(Gr+Tr)    % 遍历距离维度
    for j = Td+Gd+1 : Nd-(Gd+Td)  % 遍历多普勒维度

        noise_level = zeros(1,1);  % 初始化当前窗口的噪声能量

       
        for p = i-(Tr+Gr): i+ (Tr+Gr) % 遍历当前CUT周围的训练单元和保护单元
            for q = j-(Td+Gd): j+(Td+Gd)

                
                if (abs(i-p)> Gr ||abs(j-q)>Gd)% 排除保护单元，只计算训练单元的能量
                   % 将dB值转换为线性功率值
                    noise_level = noise_level+ db2pow(RDM(p,q));
                end
            end

        end
        % 平均训练单元能量并转换回dB
        threshold = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
        threshold = threshold + offset;
        
        % 进行检测：若CUT大于阈值，标记为目标（1），否则为噪声（0）
        CUT= RDM(i,j);
        if (CUT<threshold)
            RDM(i,j)=0;
        else 
            RDM(i,j)= 1; % max_T
            disp(i);
            disp(j);
        end
        
    end
end
% 上述过程将生成一个经过阈值处理的块，该块比距离 - 多普勒图小，
% 因为 “CUT” 不能位于矩阵的边缘。因此，有少量单元不会经过阈值处理。
% 为了保持地图大小不变，将这些单元的值设为 0。

RDM(RDM~=0 & RDM~=1) = 0;
%RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows
%RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns


% *%TODO* :
%绘制CFAR输出
figure('Name', 'CFAR')
surf(doppler_axis,range_axis,RDM);
colorbar;
title('offset 1.4');
