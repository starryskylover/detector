
PROGRAM TELESCOPEALPHA
IMPLICIT NONE
INTEGER I, J, N                        !I - number of layers; J - number of detectors
REAL E, E2, DE, DE2, DW, SW, S, POTERI
REAL W(5)   !layer thickness, мкм
REAL WSI(4) !silicon thickness, мкм
REAL A(5)   !atomic number
REAL Z(5)   !charge
REAL R(5)   !density
REAL SE(4)  !energy loss on a layer of thickness W
DATA W      /  0.01, 0.04,  99.92,  0.04,   0.02/
DATA WSI    /  9.92, 99.92, 999.92, 999.92      /
DATA A      /  197,  28,    28,     28,     27  /
DATA Z      /  79,   14,    14,     14,     13  /
DATA R      /  19.3, 2.3,   2.3,    2.3,    2.7 /
DATA SE     /  0,    0,     0,      0           /


	10 CONTINUE
	
	WRITE(*,*) 'ENTER ENERGY FROM 1 TO 40 MeV'
	READ *, E 
	
	

	
I = 1; J = 1;   
DW = W(1);      !elementary thickness in the first approximation
DE = E;         !losses in the elementary section 
SE = 0; SW = 0; !elementary thickness into which summarize
                !losses at the corresponding thickness, also summarize 

   
        DETECTOR: DO WHILE (J .LT. 5) !cycle through detectors
        
        LAYER: DO WHILE (I .LT. 6) !cycle through layers
            IF (I .EQ. 3) THEN
                W(I)=WSI(J) 
            END IF
            DW=W(I) 
            
            
            
            BLOK: DO WHILE (SW .LT. W(I)) !a separate block that takes into account 9 and 3 in algorithm
                CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)
                DE =ABS(POTERI*R(I)*DW/10000)
                
                FOURTHBLOK: DO WHILE(DE .GT. (0.05*E)) !block 4
                    DW = DW/2
                    CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)!check on the elementary thickness in the next approximation
                    DE =ABS(POTERI*R(I)*DW/10000)
                    !WRITE(*,*) 'DE=',DE
                    !WRITE(*,*)'DW=',DW
                END DO FOURTHBLOK
                
                CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)
                DE =ABS(POTERI*R(I)*DW/10000) !losses from initial energy
                E2 = E-DE   
                CALL ENERGYLOSSES(E2, Z(I), A(I), POTERI) !calculate the losses from the remaining energy 
                DE2 =ABS(POTERI*R(I)*DW/10000)
                DE = (DE + DE2)/2

                IF (I .EQ. 3) THEN
                    IF (E .LT. 0.8) THEN
                        SE(J) = SE(J) + E   !losses + residual energy
                        E = 0
                        !WRITE(*,*) 'Energy is very small, the particle has stopped'
                        EXIT DETECTOR !energy is very small, the particle has stopped
                    END IF
                    SE(J) = SE(J) + DE  !block 14
                END IF
                E = E - DE !6-8
                SW = SW + DW    !new elemental thickness
            END DO BLOK
            
            
            SW = SW - W(I)      !thickness already in the new layer 
            I = I + 1 !transition to the next layer
        END DO LAYER
        J = J +1 !move to the next detector
        I = 1 !from the first layer to the next detector
    END DO DETECTOR
    

 WRITE(*,*) 'Energy losses in detectors:'
DO N = 1, J 
    WRITE(*,*) 'E',N,'=',SE(N) 

END DO
GO TO 10
END PROGRAM TELESCOPEALPHA

SUBROUTINE ENERGYLOSSES(E, Z, A, POTERI)
    IMPLICIT NONE
    REAL E, Z, A, POTERI
    REAL CLAG1, CLAG2, CLAG3
    CLAG1=0.307*4*(Z/A)*((1+(E/3727.379))**2)/((1+(E/3727.379))**2-1)
    CLAG2=LOG(((6.2*(10000)/Z)*(E/3727.379))*((1+(0.5*E/3727.379))))
    CLAG3=1/((1+E/3727.379)**2)
    POTERI = CLAG1*(CLAG2+CLAG3)
END SUBROUTINE ENERGYLOSSES

