using PyPlot;
using Jurvis; 


function frequency(z, dt)
    y = real(z);
    ỹ = imag(z);
    ẏ = diff14(y, dt);
    ỹ̇ = diff14(ỹ, dt);

    A_sq = [(y[i]^2 + ỹ[i]^2) for i in eachindex(z)];
    #ω = Vector{Float64}(undef, length(time));
    ω = similar(y);
    for i in eachindex(ω)
        ω[i] = (y[i]*ỹ̇[i] - ẏ[i]*ỹ[i])/A_sq[i];
    end
    return ω
end



# simulation interval
dt = 1.0e-4;
t_span = [0., 5.];
t = t_span[1]:dt:t_span[2];
##
# define damping constants and time relation
m = [0.3, 0.2]; Dₐ = [0.004, 0.001]; kᵧ = [0.4, 0.8];
#damp2(x) = Dₐ[2]*(exp((x-6)/kᵧ[2])+1)-0.0016;
damp1(x) = Dₐ[1]*(-exp((-x-0.5)/kᵧ[1])+1);
damp2(x) = Dₐ[2]*(exp((-x-0.3)/kᵧ[2])+1)-0.0007;
#
ζ₁ = damp1.(t);
ζ₂ = damp2.(t);


# define frequenci constants and time relation
kf = [1.5,1.4];  f_amp = [4.5,4.5]; shift = [250 - f_amp[1], 310 - f_amp[2]];

freq1(x) = shift[1] + f_amp[1]*(1 + exp(-kf[1]*x));
#f = shift .+ f_amp.*(1 .- exp.(-kf.*t));
f₁ = freq1.(t);
ω₁ = 2π*f₁;
ϕ₁ = cumsum(ω₁)*dt;
ϕ₁2 = cumsum(ω₁.*ζ₁)*dt;
A₁=100 .* exp.(-ϕ₁2);
y₁ = A₁ .* sin.(ϕ₁);

freq2(x) = shift[2] + f_amp[2]*(1 - exp(-kf[2]*x));
#f = shift .+ f_amp.*(1 .- exp.(-kf.*t));
f₂ = freq2.(t);
ω₂ = 2π*f₂;
ϕ₂ = cumsum(ω₂)*dt;
ϕ₂2 = cumsum(ω₂.*ζ₂ )*dt;
A₂=100 .* exp.(-ϕ₂2);
y₂ = A₂ .* sin.(ϕ₂);

ζ₃ = 0.001;
f₃ = 430;
ω₃ = 2π*f₃;
#ϕ₃ = cumsum(ω₃)*dt;
#ϕ₃2 = cumsum(ω₃.*ζ₃)*dt;
A₃=100 .* exp.((-ω₃*ζ₃).*t);
y₃ = A₃ .* sin.(ω₃.*t);



# plot(y₃)
##
# plot of initial signal's properties
figure()
subplot(2,2,1);

plot(A₁, f₁);
 #ylim([0.003,0.005])
xscale("log");
ylabel("1st mode frequency")

subplot(2,2,2);
plot(A₁, ζ₁);
 #ylim([0.003,0.005])
 ylabel("1st mode damping")
xscale("log");
yscale("log");

subplot(2,2,3);
plot(A₂, f₂);
 #ylim([0.003,0.005])
 ylabel("2nd mode frequency")
xscale("log");

subplot(2,2,4);
plot(A₂, ζ₂);
 #ylim([0.003,0.005])
 ylabel("2nd mode damping")
xscale("log");
yscale("log");

tight_layout()
##

# time plots of separated modes and summarized signal

Y=y₁+y₂+y₃;
figure()
ax1 = subplot(4,1,1)
plot(y₁)
subplot(4,1,2, sharex=ax1)
plot(y₂)
subplot(4,1,3, sharex=ax1)
plot(y₃)
subplot(4,1,4,sharex=ax1)
plot(Y)

gcf().set_size_inches(6, 12)

tight_layout()

##
T=easyspectrum(Y;dt);
T1=easyspectrum(y₁;dt);
T2=easyspectrum(y₂;dt);
T3=easyspectrum(y₃;dt);
# subplot(4,1,1)
plot(T1[1],T1[2])
yscale("log")
xlim([200,500])
# subplot(4,1,2)
plot(T2[1],T2[2])
xlim([200,500])
yscale("log")
# subplot(4,1,3)
plot(T3[1],T3[2])
xlim([200,500])
yscale("log")
# subplot(4,1,4)
plot(T[1],T[2])
xlim([200,500])
yscale("log")
gcf().set_size_inches(6, 4)

tight_layout()
##
SSA_win = 6325;
signal = MeasuredData("test data", t, hcat(Y, y₁, y₂,y₃), "time", ["test signal", "mode1", "mode2","mode3",], size(Y,1), 4, ["", "", "",""]);
decomposed_data = decomposeSSA(signal,1, 20, SSA_win);

plotall(spectrumall(decomposed_data)); yscale("log"); legend();xlim([200,500]);

##

groups = [  [10,11,14, 15, 18, 19],
            [2,3,4,5,8,9],
            [6,7]
            ];
grouped_modes = groupmodes(decomposed_data, groups);#в канале 1 исходный сигнал, в каналах 2 3 4 моды 1 2 3 соответственно
#plotall(spectrumall(grouped_modes)); yscale("log"); legend();xlim([200,500]);

mode1 = copydata(grouped_modes, [2]);
# plotchannel(mode1; ch = 1); plot(t, y₁, "--")
#cut!(mode1,[0., 1.2])
env1=envelope(mode1.ydata[:,1]);
# plot(env1)
env_modes1 = decomposeSSA(env1, 10, 1000);
# plot(env_modes1)
add_data!(mode1, env_modes1[:,1]; ydata_name = "Огибающая первой моды"); #канал 2
phase1 = instphase(mode1.ydata[:,1]);
add_data!(mode1, phase1; ydata_name = "inst phase моды 1"); #канал 3
freq_modes1=decomposeSSA(findiff(mode1.ydata[:,3], dt; order = 4), 10,1000)[:,1];
ωζ1 = instdamping(freq_modes1, mode1.ydata[:,2], dt);#произведение ω и ζ
Da1=ωζ1./freq_modes1;

add_data!(mode1, freq_modes1; ydata_name = "instant frequensy of mode 1") # channel 4
add_data!(mode1, Da1; ydata_name = "instant damping of mode 1") # channel 5
# cut!(mode1, [0.1,1])
# plot( Da1);plot(ζ₁); ylim([0,0.01])# НЕ сошлось
plot(mode1.ydata[:,2], mode1.ydata[:,5]);plot(A₁,ζ₁); #ylim([0,0.01])# НЕ сошлось
yscale("log"); xscale("log");
xlim([0.1, 100]);ylim([1e-3, 1e-2])
# plotchannel(mode1; ch = 1); plot(t, y₁)
# plot(T1[1], T1[2]);plotchannel(spectrumall(mode1));yscale("log")
##

#typeof(freq_modes1)

mode2 = copydata(grouped_modes, [3]);
env2=envelope(mode2.ydata[:,1]);
env_modes2 = decomposeSSA(env2, 10, 1000);
add_data!(mode2, env_modes2[:,1]; ydata_name = "Огибающая второй моды"); #канал 2
phase2 = instphase(mode2.ydata[:,1]);
add_data!(mode2, phase2; ydata_name = "inst phase моды 2"); #канал 3
freq_modes2=decomposeSSA(findiff(mode2.ydata[:,3], dt; order = 4), 10,1000)[:,1];
ωζ2 = instdamping(freq_modes2, mode2.ydata[:,2], dt);#произведение ω и ζ
Da2=ωζ2./freq_modes2;
plot(Da2);plot(ζ₂);ylim([-0.01,0.01])# сходится



mode3 = copydata(grouped_modes, [4]);
env3=envelope(mode3.ydata[:,1]);
env_modes3 = decomposeSSA(env3, 10, 1000);
add_data!(mode3, env_modes3[:,1]; ydata_name = "Огибающая третьей моды"); #канал 2
phase3 = instphase(mode3.ydata[:,1]);
add_data!(mode3, phase3; ydata_name = "inst phase моды 3"); #канал 3
freq_modes3=decomposeSSA(findiff(mode3.ydata[:,3], dt; order = 4), 10,1000)[:,1];
ωζ3 = instdamping(freq_modes3, mode3.ydata[:,2], dt);#произведение ω и ζ
Da3=ωζ3./freq_modes3;
ζ3=ζ₃*ones(length(Da3),1);
plot(Da3);plot(ζ3);ylim([-0.01,0.01])#СОШЛОСЬ ζ₃ = 0.001;

##

#нахождение мод вычитанием из исходного сигнала с преминением определенного размера окна 
SIGNAL = MeasuredData("test data", t, hcat(Y, y₁, y₂,y₃), "time", ["test signal", "mode1", "mode2","mode3",], size(Y,1), 4, ["", "", "",""]);
plotall(spectrumall(SIGNAL)); yscale("log"); legend();xlim([200,500]);
f = [253,308,430];#по горафику определяем примерную частоту
SSA_win2 = trunc(Int,round(100/(dt*f[2])));
decomposed_data2 = decomposeSSA(SIGNAL,1, 20, SSA_win2);
plotall(spectrumall(decomposed_data)); yscale("log"); legend();xlim([200,500]);
groups2 = [  [2,3,6,7,10,11,14,15,18,19],
           
            ];
grouped_modes2 = groupmodes(decomposed_data2, groups2);#в канале 1 исходный сигнал, в каналe 2 мода 2 соответственно
plotall(spectrumall(grouped_modes2)); yscale("log"); legend();xlim([200,500]);#проверочный график моды и исходного сигнала
SminusMode2=SIGNAL.ydata[:,1]-grouped_modes2.ydata[:,2];

SIGNAL2 = MeasuredData("test data", t, hcat(SminusMode2), "time", ["signal minus mode 2" ], size(Y,1), 1, ["",]);
SSA_win3 = trunc(Int,round(100/(dt*f[3])));
decomposed_data3 = decomposeSSA(SIGNAL2,1, 10, SSA_win3);
plotall(spectrumall(decomposed_data3)); yscale("log"); legend();xlim([200,500]);
groups3 = [  [2,3]];
grouped_modes3 = groupmodes(decomposed_data3, groups3);#в канале 1 исходный сигнал, в каналe 2 мода 3 соответственно
plotall(spectrumall(grouped_modes3)); yscale("log"); legend();xlim([200,500]);#проверочный график моды и исходного сигнала
SminusMode23=SIGNAL2.ydata[:,1]-grouped_modes3.ydata[:,2];

SIGNAL3 = MeasuredData("test data", t, hcat(SminusMode23), "time", ["signal minus mode 2 3" ], size(Y,1), 1, [""]);# оставшаяся мода 1 
plotall(spectrumall(SIGNAL3)); yscale("log"); legend();xlim([200,500]);

##
Mode1 = copydata(SIGNAL3, [1]);
# plotchannel(mode1; ch = 1); plot(t, y₁, "--")
#cut!(Mode1,[0., 1.2])
env1=envelope(Mode1.ydata[:,1]);
# plot(env1)
env_modes1 = decomposeSSA(env1, 10, 1000);
# plot(env_modes1)
add_data!(Mode1, env_modes1[:,1]; ydata_name = "Огибающая первой моды"); #канал 2
phase1 = instphase(Mode1.ydata[:,1]);
add_data!(Mode1, phase1; ydata_name = "inst phase моды 1"); #канал 3
freq_modes1=decomposeSSA(findiff(Mode1.ydata[:,3], dt; order = 4), 10,1000)[:,1];
ωζ1 = instdamping(freq_modes1, Mode1.ydata[:,2], dt);#произведение ω и ζ
Da1=ωζ1./freq_modes1;

add_data!(Mode1, freq_modes1; ydata_name = "instant frequensy of mode 1"); # channel 4
add_data!(Mode1, Da1; ydata_name = "instant damping of mode 1"); # channel 5
# cut!(mode1, [0.1,1])
# plot( Da1);plot(ζ₁); ylim([0,0.01])# НЕ сошлось
plot(Mode1.ydata[:,2], Mode1.ydata[:,5]);plot(A₁,ζ₁); #ylim([0,0.01])# НЕ сошлось
yscale("log"); xscale("log");
xlim([0.1, 100]);ylim([1e-3, 1e-2])