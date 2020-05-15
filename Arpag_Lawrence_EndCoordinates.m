
% Arpag, Lawrence et al., PNAS, 2020
%
% This code is used for determining microtubule end positions and to correct for gliding of two-color kymographs.
% User first needs to open a kymograph stack in ImageJ and click on beginning and end of each growth and shrinkage events.
% The coordinates of the clicked points need to be saved as csv files with filenames in the form of 
% "Results-*-plus.csv" for plus end points,
% "Results-*-minus.csv" for minus end points.
% Similarly, user needs to click two points on kyographs to correct for
% gliding. Fiducial marks such as two-color interface can be used to
% determine the gliding. Coordinates of these points needs to be saved as
% "Results-*-gliding.csv".
% This code will then read the kymograph stack and "Results-*" files.
% Input parameters:
% "pathtofile" is foldername which contains the kymograph stack and the coordinate files.
% "prefix" corresponds to * in the filenames described above.
% "pospix" is temporal resolution.
% "timepix" is the frame interval.
% "kymoname" is the name of the kymograph stack.
% EXAMPLE USAGE:
% For a kymograph stack named as "test.tif", csv files named as
% Results-test-plus.csv, Results-test-minus.csv, Results-test-gliding.csv
% placed in the same folder as this code, run:
%%% Arpag_Lawrence_GlidingCorrection('./.','test',65,15,'test.tif') %%%
% where 65 represents 65nm temporal resolution, 15 represents 15s frame
% interval.

function Arpag_Lawrence_EndCoordinates(pathtofile,prefix,pospix,timepix,kymoname)
try
    
    close all
    clearvars -except pathtofile prefix pospix timepix kymoname
    
    newdirname=sprintf('%s/tracks/',pathtofile);
    if (0==exist(newdirname,'dir'))
        newdirname=sprintf('%s/tracks/',pathtofile);
        mkdir(newdirname)
    end
    
    kymostackname=sprintf('%s/%s',pathtofile,kymoname);
    
    %%% Read the excel file that contains plus end coordinates
    excelnameplus=sprintf('%s/Results-%s-plus.csv',pathtofile,prefix);
    datareadplus = readtable(excelnameplus);

    %%% Read the excel file that contains minus end coordinates
    excelnameminus=sprintf('%s/Results-%s-minus.csv',pathtofile,prefix);
    datareadminus = readtable(excelnameminus);
    
    %%% Read the excel file that contains gliding information
    excelnamegliding=sprintf('%s/Results-%s-gliding.csv',pathtofile,prefix);
    datareadgliding = readtable(excelnamegliding);
    
    %maxslice=datareadplus.Frame(end);
    maxslice=datareadplus.Slice(end);
    
    
    if (maxslice<datareadminus.Slice(end))
        maxslice=datareadminus.Slice(end);
    end
    if (maxslice<datareadgliding.Slice(end))
        maxslice=datareadgliding.Slice(end);
    end
    fprintf('max slice count=%d\n',maxslice)
    
    vgliding=NaN(maxslice,1);
    
    glidingcount=0;
    %%% first get the gliding information
    for stackslice=1:maxslice
        for i=1:2:datareadgliding.Var1(end)
            if (datareadgliding.Slice(i)==stackslice)
                
                
                tg1=datareadgliding.Y(i+0); % in pixel
                tg2=datareadgliding.Y(i+1); % in pixel
                
                xg1=datareadgliding.X(i+0); % in pixel
                xg2=datareadgliding.X(i+1); % in pixel
                
                glidingcount=glidingcount+1;
                vgliding(stackslice) = abs(xg2-xg1)/(tg2-tg1); % in pixel
                
                %fprintf('slice=%d\tvgliding=%f nm/s\t%d\n',stackslice,vgliding(glidingcount),glidingcount)
            end
        end
    end
    avevgliding=nanmean(vgliding); % in pixel
    stdvgliding=nanstd(vgliding); % in pixel
    fprintf('ave gliding speed: %f\t%f (nm/s)\n',avevgliding*pospix/timepix,stdvgliding*pospix/timepix)
    
    
    
    [K,map]=imread(kymostackname,1);
    
    [Nt,Nx,Nc]=size(K);
    alllength=zeros(Nt,3*maxslice); % columns will be t l lcorrected, t l lcorrected etc
    
    
    %%% now read the plus and minus end coordinates
    
    for stackslice=1:maxslice
        
        plusend=zeros(Nt,3); % time, pos, gliding-corrected-pos
        minusend=zeros(Nt,3); % time, pos, gliding-corrected-pos
        K=imread(kymostackname,stackslice);
        
        %%% determine whether plus end is on the left or right
        firstxplus=NaN;
        firstxminus=NaN;
        for i=1:datareadplus.Var1(end)-1
            if (datareadplus.Slice(i)==stackslice)
                firstxplus=datareadplus.X(i);
                break;
            end
        end
        for i=1:datareadminus.Var1(end)-1
            if (datareadminus.Slice(i)==stackslice)
                firstxminus=datareadminus.X(i);
                break;
            end
        end
        
        if ( isnan(firstxplus) || isnan(firstxminus) )
            fprintf('no plus end or minus end data found on slice %d\nSkipping..\n',stackslice)
            continue;
        end
        
        for i=1:datareadplus.Var1(end)-1
            if ( (datareadplus.Slice(i+0)==stackslice) && (datareadplus.Slice(i+1)==stackslice) )
                toi1=datareadplus.Y(i+0);
                toi2=datareadplus.Y(i+1);
                xoi1=datareadplus.X(i+0);
                xoi2=datareadplus.X(i+1);
                
                for j=1:Nt
                    if ( (j>=toi1) && (j<toi2) )
                        
                        deltat=j-toi1; % in pixels
                        v=(xoi2-xoi1)/(toi2-toi1); % in pixels
                        ep=xoi1+(v*deltat); % in pixels
                        if (~isnan(vgliding(stackslice,1)))
                            vglidingtemp=vgliding(stackslice,1);
                        else
                            vglidingtemp=avevgliding;
                        end
                        if (firstxplus>=firstxminus) % plus is on the right
                            vcorrected=1.0*vglidingtemp;
                        else
                            vcorrected=-1.0*vglidingtemp;
                        end
                        epcorrected=ep+(vcorrected*j);

                        plusend(j,1)=deltat+toi1;
                        plusend(j,2)=ep;
                        plusend(j,3)=epcorrected;
                        
                    end
                end
            end
        end
        
        
        
        for i=1:datareadminus.Var1(end)-1
            if ( (datareadminus.Slice(i+0)==stackslice) && (datareadminus.Slice(i+1)==stackslice) )
                toi1=datareadminus.Y(i+0);
                toi2=datareadminus.Y(i+1);
                xoi1=datareadminus.X(i+0);
                xoi2=datareadminus.X(i+1);
                
                for j=1:Nt
                    if ( (j>=toi1) && (j<toi2) )
                        
                        deltat=j-toi1; % in pixels
                        v=(xoi2-xoi1)/(toi2-toi1); % in pixels
                        ep=xoi1+(v*deltat); % in pixels
                        if (~isnan(vgliding(stackslice,1)))
                            vglidingtemp=vgliding(stackslice,1);
                        else
                            vglidingtemp=avevgliding;
                        end
                        if (firstxplus>=firstxminus) % plus is on the right
                            vcorrected=1.0*vglidingtemp;
                        else
                            vcorrected=-1.0*vglidingtemp;
                        end
                        
                        epcorrected=ep+(vcorrected*j);
                        
                        minusend(j,1)=deltat+toi1;
                        minusend(j,2)=ep;
                        minusend(j,3)=epcorrected;
                        
                    end
                end
            end
        end
        
        fncoord=sprintf('%s/tracks/%s-coordinates-%.3d.dat',pathtofile,prefix,stackslice);
        fidcoord=fopen(fncoord,'w');
        fncoordcorrected=sprintf('%s/tracks/%s-coordinatescorrected-%.3d.dat',pathtofile,prefix,stackslice);
        fidcoordcorrected=fopen(fncoordcorrected,'w');
        
        if (size(minusend,1)==size(plusend,1)) % this line is just to double check, no actual need
            for i=1:size(plusend,1)
                
                if ( (plusend(i,1)~=0) && (minusend(i,1)~=0) )   
                    
                    alllength(i,((3*stackslice)-2))=plusend(i,1);
                    alllength(i,((3*stackslice)-1))=abs(plusend(i,2)-minusend(i,2));
                    alllength(i,((3*stackslice)-0))=abs(plusend(i,3)-minusend(i,3));
                    
                    if (firstxplus>=firstxminus) % plus is on the right
                        
                        fprintf(fidcoord,'%f\t%f\t%f\t%f\t%f\n',timepix*plusend(i,1),pospix*minusend(i,2),pospix*plusend(i,2),pospix*alllength(i,((3*stackslice)-1)),pospix*(plusend(i,2)+minusend(i,2))/2);
                        fprintf(fidcoordcorrected,'%f\t%f\t%f\t%f\t%f\n',timepix*plusend(i,1),pospix*minusend(i,3),pospix*plusend(i,3),pospix*alllength(i,((3*stackslice)-0)),pospix*(plusend(i,3)+minusend(i,3))/2);
                    else
                        
                        fprintf(fidcoord,'%f\t%f\t%f\t%f\t%f\n',timepix*plusend(i,1),pospix*(Nx-minusend(i,2)),pospix*(Nx-plusend(i,2)),pospix*alllength(i,((3*stackslice)-1)),pospix*(plusend(i,2)+minusend(i,2))/2);
                        fprintf(fidcoordcorrected,'%f\t%f\t%f\t%f\t%f\n',timepix*plusend(i,1),pospix*(Nx-minusend(i,3)),pospix*(Nx-plusend(i,3)),pospix*alllength(i,((3*stackslice)-0)),pospix*(plusend(i,3)+minusend(i,3))/2);
                    end
                    
                else
                    
                    fprintf(fidcoord,'%f\t%f\t%f\t%f\t%f\n',0,0,0,0,0);
                    fprintf(fidcoordcorrected,'%f\t%f\t%f\t%f\t%f\n',0,0,0,0,0);
                end
                
            end
            
        end
        
        fclose(fidcoord);
        fclose(fidcoordcorrected);

        disp('ok')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% correct the images for gliding %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [K,map]=imread(kymostackname,stackslice);
        %K(:,all(all(K == 0,3),1),:) = []; % remove zero or black pixels around kymograph
        if (~isnan(vgliding(stackslice,1)))
            vglidingtemp=vgliding(stackslice,1);
        else
            vglidingtemp=avevgliding;
        end
        if (firstxplus>=firstxminus) % plus is on the right
            vcorrected=1.0*vglidingtemp;
        else
            vcorrected=-1.0*vglidingtemp;
        end
        
        figure
        subplot(1,2,1)
        imshow(K,map)
        hold on
        
        % first create
        subpixel=1;
        %Ksubpixel=zeros(Nt,subpixel*Nx,Nc);
        Kshifted=K;
        for ti=1:Nt
            for ci=1:Nc
                for xi=1:Nx
                    for xsubi=1:subpixel
                        subindex=((xi-1)*subpixel)+xsubi;
                        Kshifted(ti,subindex,ci)=0;
                        shift=floor(xi-vcorrected*ti);
                        if ( (shift>=1) && (shift<=Nx) )
                            Kshifted(ti,subindex,ci)=K(ti,shift,ci);
                        else
                            Kshifted(ti,subindex,ci)=0;
                        end
                    end
                end
            end
        end
        
        Kshifted(:,all(all(Kshifted == 0,3),1),:) = []; % remove zero or black pixels around kymograph
        subplot(1,2,2)
        imshow(Kshifted,map)
        imwrite(Kshifted,sprintf('Corrected-%d.tif',stackslice))
        
        hold on
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% FINISH: correct the images for gliding %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    
catch fnme
    fprintf('ERROR: %s\n',fnme.message)
    
end
end