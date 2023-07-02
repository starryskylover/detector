
PROGRAM TELESCOPEALPHA
IMPLICIT NONE
INTEGER I, J, N                        !I - количество слоёв; J - количество детекторов
REAL E, E2, DE, DE2, DW, SW, S, POTERI
REAL W(5)   !толщина слоев, мкм
REAL WSI(4) !толщина кремниевого ППД, мкм
REAL A(5)   !атомное число
REAL Z(5)   !заряд
REAL R(5)   !плотность
REAL SE(4)  !потери энергии на слое толщиной W
DATA W      /  0.01, 0.04,  99.92,  0.04,   0.02/
DATA WSI    /  9.92, 99.92, 999.92, 999.92      /
DATA A      /  197,  28,    28,     28,     27  /
DATA Z      /  79,   14,    14,     14,     13  /
DATA R      /  19.3, 2.3,   2.3,    2.3,    2.7 /
DATA SE     /  0,    0,     0,      0           /


	10 CONTINUE
	
	WRITE(*,*) 'ВВЕДИТЕ ЭНЕРГИЮ ОТ 1 ДО 40 МэВ'
	READ *, E 
	
	

	
I = 1; J = 1;   !I - количество слоёв; J - количество детекторов
DW = W(1);      !элементарная толщина в первом приближении
DE = E;         !потери в эл участке 
SE = 0; SW = 0; !элементарная толщина, в котрую суммируем
                !потери на соответствующей толщине, также суммируем 

   
        DETECTOR: DO WHILE (J .LT. 5) ! по числу детекторов
        
        LAYER: DO WHILE (I .LT. 6) ! по числу слоёв 
            IF (I .EQ. 3) THEN
                W(I)=WSI(J) !чтобы не вводить 20 значений для W
            END IF
            DW=W(I) !элементарная толщина слоя в первом приближении 
            
            
            
            BLOK: DO WHILE (SW .LT. W(I)) !отдельный блок, который учитывает 9 и 3
                CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)
                DE =ABS(POTERI*R(I)*DW/10000)
                
                FOURTHBLOK: DO WHILE(DE .GT. (0.05*E)) !блок 4, то есть в 9 попадют только при выполнении 4
                    DW = DW/2
                    CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)!проверяем на элементарной толщине в следующем приближении
                    DE =ABS(POTERI*R(I)*DW/10000)
                    !WRITE(*,*) 'DE=',DE
                    !WRITE(*,*)'DW=',DW
                END DO FOURTHBLOK
                
                CALL ENERGYLOSSES(E, Z(I), A(I), POTERI)
                DE =ABS(POTERI*R(I)*DW/10000) !потери от начальной энергии
                E2 = E-DE   !начальная энергия-потери
                CALL ENERGYLOSSES(E2, Z(I), A(I), POTERI) !вычисляем потери от оставшейся энергии 
                DE2 =ABS(POTERI*R(I)*DW/10000)
                DE = (DE + DE2)/2

                IF (I .EQ. 3) THEN
                    IF (E .LT. 0.8) THEN
                        SE(J) = SE(J) + E   !потери+остаточная энергия
                        E = 0
                        !WRITE(*,*) 'Оставшаяся энергия очень маленькая, частица остановилась'
                        EXIT DETECTOR !оставшаяся энергия очень маленькая, частица остановилась
                    END IF
                    SE(J) = SE(J) + DE  !блок 14
                END IF
                E = E - DE !6-8
                SW = SW + DW    !новая элементраная толщина
            END DO BLOK
            
            
            SW = SW - W(I)      ! толщина уже в новом слое 
            I = I + 1 !переход на следующий слой
        END DO LAYER
        J = J +1 !переход к следующему детектору
        I = 1 !с первого слоя для следующего детектора
    END DO DETECTOR
    

 WRITE(*,*) 'Потери энергии в детекторах:'
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

