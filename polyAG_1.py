#-*- coding: utf -*-*- 
from math import *

GroupNames = ['─CH2─','─CH(CH3)─','─CH(C6H5)─','─C(CH3)2─','─CH=CH─','─CH=C(CH3)─','─C6H4─','─O─','─CO─','─CH(OH)─','─CH(CN)─','─CO─NH─','─CF2─','─CHCl─','─CH=CCl─']
PropNames = ['Z','M','Vg','Vr','Vw','Cps','Cpl','dHm','Ek','Yg','Ym','U','Pv','δ']
Incr = []
Incr.append([1,14.03,15.85,16.45,10.23,6.05,7.26,0.90,1.00,170,170,895,20.64,0.50])
Incr.append([1,28.05,33.35,32.65,20.45,11.10,13.80,1.50,2.40,336,742,1821,41.16,0.55])
Incr.append([1,90.12,82.15,74.50,52.62,24.16,34.40,5.89,6.90,576,856,4505,147.00,0.85])
Incr.append([1,42.08,52.40,50.35,30.67,16.23,19.36,2.01,3.00,226,380,2762,61.72,0.60])
Incr.append([2,26.04,27.70,27.75,16.94,8.90,10.20,0.50,2.00,340,597.64,1490,6.32,0.75])
Incr.append([2,40.06,43.26,42.80,27.16,14.33,17.70,1.10,2.70,480,864.26,2161,7.02,0.80])
Incr.append([4,76.09,65.50,61.40,43.32,18.80,27.00,5.30,6.00,1850,2100,3556,128.60,2.50])
Incr.append([1,16,10.00,8.50,5.25,4.02,8.50,0.40,1.50,280,450,348,30.00,0.35])
Incr.append([1,28.01,13.40,19.22,11.70,5.50,12.60,1.55,2.44,365,645.25,872,65.00,0.65])
Incr.append([1,30.03,19.15,22.66,14.18,7.77,15.70,1.63,2.54,438.90,892,1088,53.50,0.70])
Incr.append([1,39.04,28.95,32.77,21.48,9.70,14.21,5.89,6.90,476.39,857.39,1768,73.50,0.90])
Incr.append([2,43.03,24.90,30.11,19.56,11.00,21.50,0.70,14.50,1000,2560,1367,125.00,1.25])
Incr.append([1,50.01,26.40,24.75,15.33,11.70,11.80,0.17,0.80,300,718,1640,70.00,1.44])
Incr.append([1,48.48,29.35,28.25,19.00,10.18,14.50,2.21,3.20,538,946,1520,90.00,0.70])
Incr.append([2,60.49,41.07,38.40,25.72,13.41,18.40,2.51,3.50,828.28,1527.52,2060,9.28,1.20])

GroupsList = []
lenGroups = len(GroupNames)

for i in range(0,lenGroups):
    GroupsList.append(0)
     
krist=0.3
Temp = 298.0
Mavg = 1.35E6
def CalcProp(GroupsList,krist,Temp,Mavg=0):
    if GroupsList[11]>0:
        Incr[0][9] = 270
    if GroupsList[6]>0 and GroupsList[11]>0:
        Incr[6][9] = 2200
        Incr[6][10] = 2300
        Incr[7][10] = 600
    
    Props = {}
    Z = 0
    M = 0
    Vg = 0
    Vr = 0
    Vw = 0
    Yg = 0
    Ym = 0
    Cps = 0
    Cpl = 0
    dHm = 0
    Ek = 0
    U = 0
    Pv = 0
    delt = 0
    for i in range(lenGroups):
        Z += Incr[i][0]*GroupsList[i]
        
    if Z<=4: Incr[0][10] = 170
    else: Incr[0][10] = 380
       
    for i in range(lenGroups):
        M += Incr[i][1]*GroupsList[i]
        Vg += Incr[i][2]*GroupsList[i]
        Vr += Incr[i][3]*GroupsList[i]
        Vw += Incr[i][4]*GroupsList[i]
        Yg += Incr[i][9]*GroupsList[i]
        Ym += Incr[i][10]*GroupsList[i]
        Cps += Incr[i][5]*GroupsList[i]
        Cpl += Incr[i][6]*GroupsList[i]
        dHm += Incr[i][7]*GroupsList[i]
        Ek += Incr[i][8]*GroupsList[i]
        U += Incr[i][11]*GroupsList[i]
        Pv += Incr[i][12]*GroupsList[i]
        delt += Incr[i][13]*GroupsList[i]
                    
    eps_g = 4.5E-4*Vw
    eps_l = 10.0E-4*Vw
    T_g = Yg/Z 
    T_m = Ym/Z 
    Vl = Vg + eps_g*(Temp-T_g) + eps_l*(Temp-T_m)
    Cp = Cps*krist + Cpl*(1-krist)
    CpsT = Cps*(0.106+0.003*Temp)
    CplT = Cpl*(0.64+0.0012*Temp)
    CpT = CpsT*krist + CplT*(1-krist)   
    cp = Cp/M
    cpT = CpT/M
    cplT = CplT/M 
    if Temp < T_m: nu = 1/3.0
    else: nu = 0.499
    
    K = (U/Vg)**6*M/Vg*1E-10
    G = 3*K*(1-2*nu)/(2*(1+nu))
    E = 3*K*(1-2*nu) 
    lam_a = 0.6E-8*cp*M/Vl*(U/Vl)**3*sqrt(3*(1-nu)/(1+nu)) 
    
    A = 27.0E3*Ek/(2*T_g**2)
    B = T_g/100 - 2
    lgn_a = 73.7*exp(-2.14*Temp/T_g) 
    lgn_k = A*(lgn_a-B)
    K_th = ((delt/Z)/sqrt(M/Z))**3
    M_cr = (0.14/K_th)**2 
    if Mavg > M_cr:
        lgn = lgn_k+3.4*log10(Mavg/M_cr)
    elif Mavg > 0:
        lgn = lgn_k-log10(M_cr/Mavg)
        
    Props['Молекулярная масса: ']=M
    Props['Плотность полимера в стеклообразном состоянии, кг/м3:  ']=M/Vg*1000.0
    Props['Плотность полимера в высокоэластическом состоянии, кг/м3:']=M/Vr*1000.0
    Props['Плотность стеклообразного полимера с учетом кристаллизации, кг/м3: ']=M/Vg*(1+0.13*krist)*1000.0
    Props['Плотность высокоэластического полимера с учетом кристаллизации, кг/м3: ']=M/Vr*(1+0.13*krist)*1000.0
    Props['Температура стеклования, К: '] = T_g
    Props['Температура плавления, К: '] = T_m
    Props['Плотность расплава, кг/м3: '] = M/Vl*1000.0
    Props['Удельный коэффициент теплового расширения в стеклообразном состоянии, м3/(кг∙К): '] = eps_g/M*1.0E-3
    Props['Удельный коэффициент теплового расширения в расплаве, м3/(кг∙К): '] = eps_l/M*1.0E-3
    Props['Удельная теплоемкость при 298 К, кДж/(кг∙К): '] = cp* 4.184
    Props['Удельная теплоемкость при заданной температуре, кДж/(кг∙К): '] = cpT* 4.184
    Props['Удельная теплоемкость расплава при заданной температуре, кДж/(кг∙К): '] = cpT* 4.184
    Props['Теплота плавления, кДж/моль: '] = dHm* 4.184
    Props['Мольная энергия когезии, кДж/моль: '] = Ek* 4.184
    
    Props['Объёмный модуль упругости при 298К, ГПа: '] = K
    Props['Модуль сдвига при заданной температуре, ГПа: '] = G
    Props['Модуль Юнга при заданной температуре, ГПа: '] = E
    
    Props['Диэлектрическая проницаемость: '] = (Pv/M)**2
    Props['Коэффициент теплопроводности, Вт/(м*К): '] = lam_a*4.184E2
    
    Props['Константа KӨ Марка-Хоувинка: '] = K_th
    Props['Критическая мольная масса полимера: '] = M_cr
    Props['Вязкость при заданной температуре, ГПа: '] = 10**lgn*1E-10

    return(Props)


import PySimpleGUI as sg

layout = []
layout.append([sg.Text('Введите количество атомных групп в звене полимера')])

for i in range(lenGroups):
    layout.append([
    sg.Text(GroupNames[i], justification='r', size=(12,1), key='-TEXT'+str(i)+'-'),
    sg.Input(size=(4,1),do_not_clear=True, key='-IN'+str(i)+'-',default_text = "0")])

layout.append([
    sg.Text("Температура, К", justification='l', size=(35,1), key='-T_Temp-'),
    sg.Input(size=(4,1),do_not_clear=True, key='-Temp-',default_text = "298")])

layout.append([
    sg.Text("Степень кристалличности, отн. ед.", justification='l', size=(35,1), key='-T_krist-'),
    sg.Input(size=(4,1),do_not_clear=True, key='-krist-',default_text = "0")])
    
layout.append([
    sg.Text("Средняя молекулярная масса", justification='l', size=(35,1), key='-T_Mavg-'),
    sg.Input(size=(12,1),do_not_clear=True, key='-Mavg-',default_text = "1.0E5")])

layout.append([sg.Button('Принять', bind_return_key=True), sg.Button('Выход')])
window = sg.Window('Расчет свойств полимеров', layout, finalize=True)

def save_results(props,GroupsList):

    with open("results.txt", 'a') as f:
        line = 'Группы:\n'
        for i in range(len(GroupsList)):
            if GroupsList[i]>0:
                line += GroupNames[i]+'(*'+str(GroupsList[i])+')'+'\n'                
        for key,value in props.items():
            if value > 1 and value < 1000:
                line += key + ' ' + str(round(value,2)) + '\n'
            else:    
                line += key + ' ' + str('{:.3e}'.format(value)) + '\n'
        line += '\n'
        f.write(line)

    sg.popup('Результаты сохранены в файле results.txt')


def props_window(props,GroupsList):
    layout = [
              [sg.Output(size=(110,24))],
              [sg.Button('Сохранить'),sg.Button('Выход')],
              ]
    window = sg.Window('Результат расчета', layout, finalize=True)
    for key,value in props.items():
        if value > 1 and value < 1000:
            print(key, round(value,2))
        else:    
            print(key, '{:.3e}'.format(value))
 
    while True:
        event, values = window.read()
        if event == 'Сохранить':
            save_results(props,GroupsList)
        if event in (None, 'Exit', 'Выход'):
           break 
    window.close()
     


while True:             # Event Loop
    event, values = window.read()
    if event in (None, 'Exit', 'Выход'):
        break
    if event == 'Принять':
        for i in range(lenGroups):
            val=values['-IN'+str(i)+'-']
            GroupsList[i]=float(val)
        Temp = float(values['-Temp-'])
        krist = float(values['-krist-'])
        Mavg = float(values['-Mavg-'])
        proplist=CalcProp(GroupsList,krist,Temp,Mavg)
        props_window(proplist,GroupsList)
       
window.close()














        




