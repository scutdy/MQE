function  MQE10
 
         num=1;
      
          n=4 ;%number of groups
          m=5;%number of particles in each group
          le=m;
          q=m*n;
  
       
   totalNUM=[1 ];
    
 for problemIndex =totalNUM
       
        
      switch problemIndex

        case 1%   sphere function
             Dmin= -100;
             Dmax =100;
           func_num=1;%评估函数
           bestvalue=0;
            M=30;%维度
        
                
      end
 
      M=30;
        MAX_FES=M*1e4;
        onebest=zeros(1,num);%每个函数25个最优解
        allbest=zeros(num,MAX_FES);%每个函授25次运行的所有优解
   
 for tdy=1:num
       
  FES=0;
  tic; t1=clock; 
  rand('seed', sum(100*clock));   
  particle=struct('fitness',{},'center',{});
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for i=1:q   
data=Dmin+(Dmax-Dmin).*rand(1,M);    
fitness=benchmark_func(data,func_num);
particle(i).fitness=fitness;
particle(i).center=data;
FES=FES+1;
 if i==1
       tempparticle=particle(1);
               
   end
           
           if particle(i).fitness< tempparticle.fitness
               tempparticle=particle(i);
             
           end
              allbest(tdy,FES)= tempparticle.fitness;
               
              fprintf('fitMQE(%d) =%g\n',FES,allbest(tdy,FES));

end
 

%%%%%%%%%%%%%%%%%
 
while FES<=MAX_FES%混合排序
        if  FES>MAX_FES
                break;
        end
  
for i=1:q-1
    for j=1:q-i
        if particle(j).fitness < particle(j+1).fitness
            temp=particle(j+1).fitness;
            temp2=particle(j+1).center;
            
             particle(j+1).fitness=particle(j).fitness;
             particle(j+1).center=particle(j).center;
             
             particle(j).center=temp2;
             particle(j).fitness=temp;
        end 
    end
end
  
  for ttt=1:m
        memory(ttt).center =zeros(size( particle(1).center,1),size( particle(1).center,2));
  end 
 

  
for i=1:n
    for k=1:le      
     %%%%%%%%%%%%%
    Xw= particle(i);
     XwNo=i;
      locXwNO=1;
      
     Xb= particle(i);
     XbNo=i;
     locXbNO=1;
     
     Xwb=particle(i);
      locMean=[];
      fitMeans=zeros(1,n);
               
            for tt=1:m
                  if  particle(i+n*(tt-1)).fitness<Xb.fitness
                      Xb=particle(i+n*(tt-1));
                       XbNo=i+n*(tt-1) ;
                      locXbNO=tt ;
                  end
                  if  particle(i+n*(tt-1)).fitness>Xw.fitness
                      Xw=particle(i+n*(tt-1));
                     XwNo=i+n*(tt-1) ;
                     locXwNO=tt;
                  end
                fitMeans(tt)=particle(i+n*(tt-1)).fitness;  
                 locMean=[locMean;
                               particle(i+n*(tt-1)).center];
            end
            
            w11=rand ;
             v11=rand ;
             Xwb.center=(w11/(w11+v11)).*Xb.center+(v11/(w11+v11)).*Xw.center;       
           
             for tt2=1:m     
                 dis(tt2)= norm(particle((i+n*(tt2-1))).center-particle(XbNo).center,2)^2 +eps;  
             end         
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
             [Mass]=massCalculation(fitMeans,1);
             
              tempDis= dis;
              tempDis([locXbNO  ])=[];
              
              tempMass= Mass;
              tempMass([locXbNO  ])=[];
                
              Gm= tempMass./tempDis;
              totalGm=sum(Gm);
              Gm=Gm./totalGm;
         
               tempLocalMeans= locMean;
               tempLocalMeans([locXbNO  ],:)=[];
               
              mbest1=  sum(locMean,1)/(n) ;
              
 
                 
              mbest2=  sum(Gm'*ones(1,M).*tempLocalMeans,1) ;
           
                if rand<0.5
                   mbest=  mbest1;
          
                else
                    mbest= mbest2;
                    
                end
                
               temp=Xw.center;
               z=rand;     
               rrr=rand(1,size(temp,2));
               NN=1.5-0.5*FES/MAX_FES;
                            
               temp1=  Xwb.center-NN*abs(mbest -particle((XwNo)).center).*log(1./z) ;
               temp2 = Xwb.center+NN*abs(mbest-particle((XwNo)).center).*log(1./z) ;
                     
               temp(find(rrr>0.5))=temp1(find(rrr>0.5)) ;        
               temp(find(rrr<=0.5))=temp2(find(rrr<=0.5))  ; 
                          
                temp(temp>Dmax)=Dmax;
                temp(temp<Dmin)=Dmin;
                           
                             
                 fitness =benchmark_func(temp,func_num);
                 FES=FES+1;
                             
                  if  fitness <particle((XwNo)).fitness
                                 particle((XwNo)).center  = temp ;
                                 particle((XwNo)).fitness=fitness;
                 end
             %记录每个FES最佳适应度值
                   if  fitness< tempparticle.fitness    
                             tempparticle.center=temp;
                             tempparticle.fitness=fitness;
                   end                               
                 if  FES>MAX_FES
                            break;
                 end
                 allbest(tdy,FES)= tempparticle.fitness;
                              
                if mod(FES,11113)==0
                               fprintf('fitMQE(%d) =%g\n',FES,allbest(tdy,FES));
                end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end%m
 
end%le



end%while
  
 toc;etime(clock,t1) ;
  onebest(tdy)=allbest (tdy,MAX_FES);

 end%tdy
 
   
 
end%21function




% sbest(tdy)
function [M]=massCalculation(fit,min_flag)
%%%%here, make your own function of 'mass calculation'

Fmax=max(fit); Fmin=min(fit); Fmean=mean(fit); 
[i N]=size(fit);

if Fmax==Fmin
   M=ones(1,N);
else
    
   if min_flag==1 %for minimization
      best=Fmin;worst=Fmax; % 
   else %for maximization
      best=Fmax;worst=Fmin; % 
   end
  
   M=1*(fit-worst)./(best-worst)+0  ; % 

end

M=M./sum(M); %eq. 16.



function f=benchmark_func(data,func_num)

switch func_num

        case 1
             
            f=f1(data);
        
end
 

   function y=f1(x)
% This is sphere function
% x is a vector[-100,100],opt=0
d=length(x);
y=0;
for k=1:d
    y=y+x(k)^2;
end
 