#### Bibliotecas ####

import math
import cmath

#### Variaveis ####

A_pol = [0,0]
C_pol = [0,0]

#### Inputs ####

print("\n-------------------------------------------------------------")
print("\n\t\tLinhas de transmissao")
print("\n-------------------------------------------------------------\n")

vn = float(input("Digite a tensao nominal de operacao da linha [kV]: "))
l = float(input("Digite o comprimento da linha [km] : "))
f = float(input("Digite a frequencia nominal da linha [Hz]: "))
w = 2*math.pi*f #[rad/s]

aux = 0

while aux!=1 and aux!=2:
	aux = int(input("\nVoce quer inserir A e C na forma polar (1) ou na retangular (2)?"))
	if aux==1:
		A_pol[0] = float(input("Digite o valor do modulo de A : "))
		A_pol[1] = float(input("Digite o valor do angulo de A : "))
		C_pol[0] = float(input("Digite o valor do modulo de C [S] : "))
		C_pol[1] = float(input("Digite o valor do angulo de C : "))

		A_rec = complex(math.cos(A_pol[1]*(math.pi/180))*A_pol[0],math.sin(A_pol[1]*(math.pi/180))*A_pol[0])
		C_rec = complex(math.cos(C_pol[1]*(math.pi/180))*C_pol[0],math.sin(C_pol[1]*(math.pi/180))*C_pol[0])

	if aux==2:
		A_rec = complex(input("Digite o valor de A na forma retangular (X+Yj) [1/km] :"))
		raio = cmath.polar(A_rec)
		angulo = raio[1]*(180/math.pi)
		A_pol = [raio[0], angulo]

		C_rec = complex(input("Digite o valor de C na forma retangular (X+Yj) [ohm] :"))
		raio = cmath.polar(C_rec)
		angulo = raio[1]*(180/math.pi)
		C_pol = [raio[0], angulo]

vrnl = float(input("\nDigite o valor da tensao nominal da carga [kV] (Vr) (0 se nao existir) : "))
p3r = float(input("Digite o valor da carga trifasica [MW] (P3r) (0 se nao existir): "))
fpr = float(input("Digite o fator de potencia FPr (0 se nao existir) : "))

#### Modelo exponencial ####

print("\n-------------------------------------------------------------")
print("\n\t\t\tModelo exponencial")
print("\n-------------------------------------------------------------\n")

gama_rec = cmath.acosh(A_rec)/l
raio = cmath.polar(gama_rec)
angulo = raio[1]*(180/math.pi)
gama_pol = [raio[0], angulo]
print("\n* O coeficiente de propagacao gama = acosh(A)/l [1/km]\n//rec:", gama_rec, "// pol:", gama_pol,"//")

gl_rec = cmath.acosh(A_rec)
raio = cmath.polar(gl_rec)
angulo = raio[1]*(180/math.pi)
gl_pol = [raio[0], angulo]
print("\n* O resultado de gama*l = acosh(A)\n//rec:", gl_rec, "// pol:", gl_pol,"//")

alfa = gama_rec.real
print("\n* O coeficiente de atenuacao alfa = re(gama) [neper/km] eh:", alfa)

beta = gama_rec.imag
print("\n* O coeficiente de fase beta = imag(gama) [rad/km] eh:", beta)

lamb = 2*math.pi/beta
print("\n* O comprimento de onda lambda = 2pi/beta [km] eh:", lamb)

v = lamb*f
print("\n* A velocidade de propagcao v = lamb*f [km/s] eh:", v)

t = l/v
print("\n* O tempo de transito t = l/v [s] eh:", t)

print("\n-------------------------------------------------------------")

zc_rec = cmath.sinh(gl_rec)/C_rec
raio = cmath.polar(zc_rec)
angulo = raio[1]*(180/math.pi)
zc_pol = [raio[0], angulo]
print("\n* A impedancia caracteristica Zc = sinh(gama*l)/C [ohm]\n// rec:", zc_rec, "// pol:", zc_pol,"//")

yc_rec = 1/zc_rec
raio = cmath.polar(yc_rec)
angulo = raio[1]*(180/math.pi)
yc_pol = [raio[0], angulo]
print("\n* A admitancia caracteristica Yc = 1/Zc [S]\n// rec:", yc_rec, "// pol:", yc_pol,"//")

#### Rs, Ls, Cs, Zs e Ys ####

print("\n-------------------------------------------------------------")
print("\n\t\t\tRs, Ls, Cs, Zs e Ys")
print("\n-------------------------------------------------------------\n")

z_rec = gama_rec/yc_rec
raio = cmath.polar(z_rec)
angulo = raio[1]*(180/math.pi)
z_pol = [raio[0], angulo]

y_rec = z_rec/cmath.sqrt(zc_rec)
raio = cmath.polar(y_rec)
angulo = raio[1]*(180/math.pi)
y_pol = [raio[0], angulo]

rs = z_rec.real
ls = z_rec.imag/w
cs = y_rec.imag/w

print("\n* A impedancia longituginal Zs = Rs+jwLs [ohm/km]\n// rec:", z_rec, "// pol:", z_pol,"//")
print("\n* A admiancia transversal Ys = jwCs [S/km]\n// rec:", y_rec, "// pol:", y_pol,"//")
print("\nO valor de Rs [ohm/km] eh:", rs)
print("\nO valor de Ls [H/km] eh:", ls)
print("\nO valor de Cs [F/km] eh:", cs)


#### SIL ####

print("\n-------------------------------------------------------------")
print("\n\t\t\tSIL")
print("\n-------------------------------------------------------------\n")

z0 = math.sqrt(ls/cs)
print("\n* A impedancia de surto Z0 [ohm] = sqrt(Ls/Cs) eh:",z0)

sil3 = (vn*vn)/z0
print("\n* O SIL trifasico [MW] = Vn^2/z0 eh:",sil3)

sil1 = sil3/3
print("\n* O SIL monofasico [MW] = sil3/3 eh:",sil1)

#### Modelo ABCD EFGH ####

print("\n-------------------------------------------------------------")
print("\n\t\t\tModelo ABCD EFGH")
print("\n-------------------------------------------------------------\n")

B_rec = zc_rec*cmath.sinh(gama_rec*l)
raio = cmath.polar(B_rec)
angulo = raio[1]*(180/math.pi)
B_pol = [raio[0], angulo]

D_rec = cmath.cosh(gama_rec*l)
raio = cmath.polar(D_rec)
angulo = raio[1]*(180/math.pi)
D_pol = [raio[0], angulo]

E_rec = A_rec
raio = cmath.polar(E_rec)
angulo = raio[1]*(180/math.pi)
E_pol = [raio[0], angulo]

F_rec = -B_rec
raio = cmath.polar(F_rec)
angulo = raio[1]*(180/math.pi)
F_pol = [raio[0], angulo]

G_rec = -C_rec
raio = cmath.polar(G_rec)
angulo = raio[1]*(180/math.pi)
G_pol = [raio[0], angulo]

H_rec =	A_rec
raio = cmath.polar(H_rec)
angulo = raio[1]*(180/math.pi)
H_pol = [raio[0], angulo]

print("\nD = cosh(gama*l) // B = Zc*sinh(gama*l) [ohm]")
print("\nE = H = A // F = -B [ohm] // G = -C [S]")
print("\nA rec:", A_rec, "\tA pol:", A_pol, "\nB rec:", B_rec, "\tB pol:", B_pol, "\nC rec:", C_rec, "\tC pol:", C_pol,"\nD rec:", D_rec, "\tD pol:", D_pol)
print("\nE rec:", E_rec, "\tE pol:", E_pol, "\nF rec:", F_rec, "\tF pol:", F_pol, "\nG rec:", G_rec, "\tG pol:", G_pol, "\nH rec:", H_rec, "\tH pol:", H_pol)

#### Modelo PI ####

print("\n-------------------------------------------------------------")
print("\n\t\t\tModelo PI")
print("\n-------------------------------------------------------------\n")

print("\n* Circuito pi equivalente para linhas longas linha>240km")
zpi_rec = B_rec
raio = cmath.polar(zpi_rec)
angulo = raio[1]*(180/math.pi)
zpi_pol = [raio[0], angulo]

ypi_rec = (A_rec-1)/B_rec
raio = cmath.polar(ypi_rec)
angulo = raio[1]*(180/math.pi)
ypi_pol = [raio[0], angulo]
print("\nO resultado Zpi [ohm] = B\n//rec:", zpi_rec,"// pol:",zpi_pol,"//")
print("\nO resultado Ypi/2 [S] = (A-1)/B\n//rec:", ypi_rec,"// pol:",ypi_pol,"//")

print("\n-------------------------------------------------------------")

print("\n* Circuito pi nominal para linhas medias 240km>linha>80km")
zn_rec=z_rec*l
raio = cmath.polar(zn_rec)
angulo = raio[1]*(180/math.pi)
zn_pol = [raio[0], angulo]

yn_rec = y_rec*l/2
raio = cmath.polar(yn_rec)
angulo = raio[1]*(180/math.pi)
yn_pol = [raio[0], angulo]
print("\nO resultado Zn [ohm] = Zs*l\n//rec:", zn_rec,"// pol:",zn_pol,"//")
print("\nO resultado Yn/2 [S] = (Ys*l)/2\n//rec:", yn_rec,"// pol:",yn_pol,"//")


print("\n-------------------------------------------------------------")

print("\n* Impedancia RL serie para linhas curtas linha<80km")
zn2_rec=z_rec*l
raio = cmath.polar(zn_rec)
angulo = raio[1]*(180/math.pi)
zn2_pol = [raio[0], angulo]

print("\nZn [ohm] = Zs*l\n//rec:", zn2_rec,"// pol:",zn2_pol,"//")

print("\n-------------------------------------------------------------")

print("\n* Relacoes")

z_rel_pol = [zpi_pol[0]/zn_pol[0],zpi_pol[1]-zn_pol[1]]

y_rel_pol = [ypi_pol[0]/yn_pol[0],ypi_pol[1]-yn_pol[1]]

z2_rel_pol = [zn_pol[0]/zpi_pol[0],zn_pol[1]-zpi_pol[1]]

y2_rel_pol = [yn_pol[0]/ypi_pol[0],yn_pol[1]-ypi_pol[1]]

print("\nO valor de Zpi/Zn pol eh:", z_rel_pol)
print("O valor de Ypi/Yn pol eh:", y_rel_pol)
print("\nO valor de Zn/Zpi pol eh:", z2_rel_pol)
print("O valor de Yn/Ypi pol eh:", y2_rel_pol)

#### Dados RECEPTOR e EMISSOR ####

if vrnl!=0 and p3r!=0:

	print("\n-------------------------------------------------------------")
	print("\n\t\t\tDados RECEPTOR de fase")
	print("\n-------------------------------------------------------------\n")

	vr_rec = complex((vrnl/cmath.sqrt(3)),0)
	raio = cmath.polar(vr_rec)
	angulo = raio[1]*(180/math.pi)
	vr_pol = [raio[0], angulo]
	print("\n* A tensão na barra receptora Vr = Vrnl/sqrt(3) [kV] na forma polar eh:", vr_pol)

	ir_modulo = p3r/(math.sqrt(3)*vrnl*fpr)
	ir_angulo = math.acos(fpr)*(180/math.pi)
	ir_pol = [ir_modulo, -ir_angulo]
	ir_rec = cmath.rect(ir_pol[0],math.radians(ir_pol[1]))
	print("\n* A corrente na barra receptora |Ir| = Pr(3phi)/(sqrt(3)*Vr*fp) com angulo phi_ir = arccos(fp) [kA] na forma polar eh:\n",ir_pol)

	sr_pol = [vr_pol[0]*ir_pol[0], vr_pol[1]+(-ir_pol[1])]
	sr_rec = cmath.rect(sr_pol[0],math.radians(sr_pol[1]))
	pr = sr_rec.real
	qr = sr_rec.imag
	print("\n* A potencia complexa na barra receptora Sr = Vr*Ir* [MVA]")
	print("Sr rec:", sr_rec, "\tSr pol:", sr_pol)
	print("Pr [MW]:", pr, "\tQr [MVAr]:", qr)

	fpr = pr/sr_pol[0]
	print("\n* O FPr eh:", fpr)

	print("\n-------------------------------------------------------------")
	print("\n\t\tDados EMISSOR de fase")
	print("\n-------------------------------------------------------------\n")

	vs_rec = (A_rec*vr_rec)+(B_rec*ir_rec)
	raio = cmath.polar(vs_rec)
	angulo = raio[1]*(180/math.pi)
	vs_pol = [raio[0], angulo]
	print("\n* A tensão na barra emissora Vs = A*Vr+B*Ir [kV] na forma polar eh:\n", vs_pol)

	is_rec = (C_rec*vr_rec)+(D_rec*ir_rec)
	raio = cmath.polar(is_rec)
	angulo = raio[1]*(180/math.pi)
	is_pol = [raio[0], angulo]
	print("\n* A corrente na barra emissora Is = C*Vr+D*Ir [kA] na forma polar eh:\n", is_pol)

	ss_pol = [vs_pol[0]*is_pol[0], vs_pol[1]+(-is_pol[1])]
	ss_rec = complex(math.cos(ss_pol[1]*(math.pi/180))*ss_pol[0],math.sin(ss_pol[1]*(math.pi/180))*ss_pol[0])
	ps = ss_rec.real
	qs = ss_rec.imag
	print("\n* A potencia complexa na barra emissora Ss = Vs*Is* [MVA]")
	print("Ss rec:", ss_rec, "\tSs pol:", ss_pol)
	print("Ps [MW]:", ps, "\tQs [MVAr]:", qs)

	fps = ps/ss_pol[0]
	print("\n* O FPs eh:", fps)

#### Indicies de desempenho ####

	print("\n-------------------------------------------------------------")
	print("\n\t\tIndicies de desempenho")
	print("\n-------------------------------------------------------------\n")

	reg = 100*((vs_pol[0]/A_pol[0])-vr_pol[0])/(vn/math.sqrt(3))
	print("\n* A regulacao de tensao Reg = 100* (|Vs/A|-|Vr(fl)|)/(Vn/sqrt(3)) eh:", reg, "%")

	N = 100*pr/ps
	print("\n* O rendimento da transmissao N = 100* Pr/Ps eh:", N, "%")

	delta_v = 100*(vs_pol[0]-vr_pol[0])/(vn/math.sqrt(3))
	print("\n* A variacao percentual de tensao delta_V = 100* (|Vs|-|Vr|)/(Vn/sqrt(3)) eh:", delta_v, "%")

	delta_q = 100*(qs-qr)/pr
	print("\n* O consumo percentual de reativo delta_Q = 100* (Qs-Qr)/Pr eh:", delta_q, "%")

#### Modelo An Bn Cn Dn ####

	print("\n-------------------------------------------------------------\n")

	aux = int(input("Voce quer o modelo An Bn Cn Dn da linha? Nao(0) Sim(1)"))

	print("\n-------------------------------------------------------------\n")

	if aux==1:

		print("\n\n-------------------------------------------------------------")
		print("\n-------------------------------------------------------------")
		print("\n\t\tModelo An Bn Cn Dn (nominal)")
		print("\n-------------------------------------------------------------")
		print("\n-------------------------------------------------------------\n")

		An_rec = 1+(zn_rec*yn_rec)
		raio = cmath.polar(An_rec)
		angulo = raio[1]*(180/math.pi)
		An_pol = [raio[0], angulo]

		Bn_rec = zn_rec
		raio = cmath.polar(Bn_rec)
		angulo = raio[1]*(180/math.pi)
		Bn_pol = [raio[0], angulo]

		Cn_rec = yn_rec+yn_rec*zn_rec*yn_rec+yn_rec
		raio = cmath.polar(Cn_rec)
		angulo = raio[1]*(180/math.pi)
		Cn_pol = [raio[0], angulo]

		Dn_rec = 1+(zn_rec*yn_rec)
		raio = cmath.polar(Dn_rec)
		angulo = raio[1]*(180/math.pi)
		Dn_pol = [raio[0], angulo]

		print("\nA = D = Zn*(Yn/2) // B = Zn [ohm] // C = (Yn/2)+(Yn/2)*Zn*(Yn/2)+(Yn/2) [S]")
		print("\nA rec:", An_rec, "\tA pol:", An_pol, "\nB rec:", Bn_rec, "\tB pol:", Bn_pol, "\nC rec:", Cn_rec, "\tC pol:", Cn_pol,"\nD rec:", Dn_rec, "\tD pol:", Dn_pol)

		print("\n-------------------------------------------------------------")
		print("\n\t\tDados RECEPTOR de fase linhas medias")
		print("\n-------------------------------------------------------------")

		print("\n* A tensão na barra receptora Vr = Vrnl/sqrt(3) [kV] na forma polar eh:", vr_pol)
		print("\n* A corrente na barra receptora |Ir| = Pr(3phi)/(sqrt(3)*Vr*fp) com angulo phi_ir = arccos(fp) [kA] na forma polar eh:\n",ir_pol)
		print("\n* A potencia complexa na barra receptora Sr = Vr*Ir* [MVA]")
		print("Sr rec:", sr_rec, "\tSr pol:", sr_pol)
		print("Pr [MW]:", pr, "\tQr [MVAr]:", qr)
		print("\n* O FPr eh:", fpr)

		print("\n-------------------------------------------------------------")
		print("\n\t\tDados EMISSOR de fase linhas medias")
		print("\n-------------------------------------------------------------")

		vsn_rec = (An_rec*vr_rec)+(Bn_rec*ir_rec)
		raio = cmath.polar(vsn_rec)
		angulo = raio[1]*(180/math.pi)
		vsn_pol = [raio[0], angulo]
		print("\n* A tensão na barra emissora Vs = An*Vr+Bn*Ir [kV] para o modelo n na forma polar eh:\n", vsn_pol)

		isn_rec = (Cn_rec*vr_rec)+(Dn_rec*ir_rec)
		raio = cmath.polar(isn_rec)
		angulo = raio[1]*(180/math.pi)
		isn_pol = [raio[0], angulo]
		print("\n* A corrente na barra emissora Is = Cn*Vr+Dn*Ir [kA] para o modelo n na forma polar eh:\n", isn_pol)

		ssn_pol = [vsn_pol[0]*isn_pol[0], vsn_pol[1]+(-isn_pol[1])]
		ssn_rec = complex(math.cos(ssn_pol[1]*(math.pi/180))*ssn_pol[0],math.sin(ssn_pol[1]*(math.pi/180))*ssn_pol[0])
		psn = ssn_rec.real
		qsn = ssn_rec.imag
		print("\n* A potencia complexa na barra emissora Ss = Vs*Is* [MVA] para o modelo n")
		print("Ss rec:", ssn_rec, "\tSs pol:", ssn_pol)
		print("Ps [MW]:", psn, "\tQs [MVAr]:", qsn)

		fpsn = psn/ssn_pol[0]
		print("\n* O FPs eh:", fpsn)

		print("\n-------------------------------------------------------------")
		print("\n\t\tIndicies de desempenho linhas medias")
		print("\n-------------------------------------------------------------")

		regn = 100*((vsn_pol[0]/An_pol[0])-vr_pol[0])/(vn/math.sqrt(3))
		print("\n* A regulacao de tensao Reg = 100* (|Vs/A|-|Vr(fl)|)/(Vn/sqrt(3)) eh:", regn, "%")

		Nn = 100*pr/psn
		print("\n* O rendimento da transmissao N = 100* Pr/Ps eh:", N, "%")

		delta_vn = 100*(vsn_pol[0]-vr_pol[0])/(vn/math.sqrt(3))
		print("\n* A variacao percentual de tensao delta_V = 100* (|Vs|-|Vr|)/(Vn/sqrt(3)) eh:", delta_vn, "%")

		delta_qn = 100*(qsn-qr)/pr
		print("\n* O consumo percentual de reativo delta_Q = 100* (Qs-Qr)/Pr eh:", delta_qn, "%")

#### Modelo de linha curta ####

	print("\n-------------------------------------------------------------\n")

	aux = int(input("Voce quer o modelo de linha curta? Nao(0) Sim(1)"))

	print("\n-------------------------------------------------------------\n")

	if aux==1:

		print("\n\n-------------------------------------------------------------")
		print("\n-------------------------------------------------------------")
		print("\n\t\tModelo de linha curta Ac Bc Dc Cc")
		print("\n-------------------------------------------------------------")
		print("\n-------------------------------------------------------------\n")

		Ac_rec = 1+0j
		raio = cmath.polar(Ac_rec)
		angulo = raio[1]*(180/math.pi)
		Ac_pol = [raio[0], angulo]

		Bc_rec = zn2_rec
		raio = cmath.polar(Bc_rec)
		angulo = raio[1]*(180/math.pi)
		Bc_pol = [raio[0], angulo]

		Cc_rec = 0+0j
		raio = cmath.polar(Cc_rec)
		angulo = raio[1]*(180/math.pi)
		Cc_pol = [raio[0], angulo]

		Dc_rec = 1+0j
		raio = cmath.polar(Dc_rec)
		angulo = raio[1]*(180/math.pi)
		Dc_pol = [raio[0], angulo]

		print("\nA = D = 1 // B = Zn [ohm] // C = 0 [S]")
		print("\nA rec:", Ac_rec, "\tA pol:", Ac_pol, "\nB rec:", Bc_rec, "\tB pol:", Bc_pol, "\nC rec:", Cc_rec, "\tC pol:", Cc_pol,"\nD rec:", Dc_rec, "\tD pol:", Dc_pol)

		print("\n-------------------------------------------------------------")
		print("\n\t\tDados RECEPTOR de fase linhas medias")
		print("\n-------------------------------------------------------------")

		print("\n* A tensão na barra receptora Vr = Vrnl/sqrt(3) [kV] na forma polar eh:", vr_pol)
		print("\n* A corrente na barra receptora |Ir| = Pr(3phi)/(sqrt(3)*Vr*fp) com angulo phi_ir = arccos(fp) [kA] na forma polar eh:\n",ir_pol)
		print("\n* A potencia complexa na barra receptora Sr = Vr*Ir* [MVA]")
		print("Sr rec:", sr_rec, "\tSr pol:", sr_pol)
		print("Pr [MW]:", pr, "\tQr [MVAr]:", qr)
		print("\n* O FPr eh:", fpr)

		print("\n-------------------------------------------------------------")
		print("\n\t\tDados EMISSOR de fase linhas medias")
		print("\n-------------------------------------------------------------")

		vsc_rec = (Ac_rec*vr_rec)+(Bc_rec*ir_rec)
		raio = cmath.polar(vsc_rec)
		angulo = raio[1]*(180/math.pi)
		vsc_pol = [raio[0], angulo]
		print("\n* A tensão na barra emissora Vs = Ac*Vr+Bc*Ir [kV] para o modelo n na forma polar eh:\n", vsc_pol)

		isc_rec = (Cc_rec*vr_rec)+(Dc_rec*ir_rec)
		raio = cmath.polar(isc_rec)
		angulo = raio[1]*(180/math.pi)
		isc_pol = [raio[0], angulo]
		print("\n* A corrente na barra emissora Is = Cc*Vr+Dc*Ir [kA] para o modelo n na forma polar eh:\n", isc_pol)

		ssc_pol = [vsc_pol[0]*isc_pol[0], vsc_pol[1]+(-isc_pol[1])]
		ssc_rec = complex(math.cos(ssc_pol[1]*(math.pi/180))*ssc_pol[0],math.sin(ssc_pol[1]*(math.pi/180))*ssc_pol[0])
		psc = ssc_rec.real
		qsc = ssc_rec.imag
		print("\n* A potencia complexa na barra emissora Ss = Vs*Is* [MVA] para o modelo n")
		print("Ss rec:", ssc_rec, "\tSs pol:", ssc_pol)
		print("Ps [MW]:", psc, "\tQs [MVAr]:", qsc)

		fpsc = psc/ssc_pol[0]
		print("\n* O FPs eh:", fpsc)

		print("\n-------------------------------------------------------------")
		print("\n\t\tIndicies de desempenho linhas medias")
		print("\n-------------------------------------------------------------")

		regc = 100*((vsc_pol[0]/Ac_pol[0])-vr_pol[0])/(vn/math.sqrt(3))
		print("\n* A regulacao de tensao Reg = 100* (|Vs/A|-|Vr(fl)|)/(Vn/sqrt(3)) eh:", regc, "%")

		Nn = 100*pr/psc
		print("\n* O rendimento da transmissao N = 100* Pr/Ps eh:", N, "%")

		delta_vc = 100*(vsc_pol[0]-vr_pol[0])/(vn/math.sqrt(3))
		print("\n* A variacao percentual de tensao delta_V = 100* (|Vs|-|Vr|)/(Vn/sqrt(3)) eh:", delta_vc, "%")

		delta_qc = 100*(qsc-qr)/pr
		print("\n* O consumo percentual de reativo delta_Q = 100* (Qs-Qr)/Pr eh:", delta_qc, "%")

print("\n-------------------------------------------------------------")
print("\n\t\t\tFim")
print("\n-------------------------------------------------------------")







