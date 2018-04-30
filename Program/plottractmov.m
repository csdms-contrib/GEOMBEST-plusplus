function plottractmov(filethread,tf,modelrun)

%function to create a movie of simulations from plottractcolour plots

% If the movie has already been made, play it
if exist(['../Output' num2str(filethread) '/pmovie.mat'])
    load(['../Output' num2str(filethread) '/pmovie.mat'])
    plottractcolour(filethread,tf,1,tf,1,modelrun)
    movie(pmovie,1,4)
else
    % Record the frames for the movies
    for t = 1:tf
        plottractcolour(filethread,t,1,t,1,modelrun)
        text(0,2.3,['t = ' num2str((t-1)*10) ' yr'],'fontsize',24)
        pmovie(t) = getframe;
        outputfilename = ['../Output' num2str(filethread) '/movieframe' num2str(t)];
        set(gcf, 'InvertHardCopy', 'off');
        print('-dpng',outputfilename)
    end
    % Play movie
    movie(pmovie,1,4)
    save(['../Output' num2str(filethread) '/pmovie.mat'],'pmovie')
end

writerObj = VideoWriter(['../Output' num2str(filethread) '/pmovie.avi']);
open(writerObj);

plottractcolour(filethread,1,1,1,1,modelrun); 
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

for t = 1:tf 
   plottractcolour(filethread,t,1,t,1,modelrun)
   text(0,2.3,['t = ' num2str((t-1)*10) ' yr'],'fontsize',24)
   frame = getframe;
   writeVideo(writerObj,frame);
end

close(writerObj);