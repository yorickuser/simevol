# file simevol/R/simevol1.R
# copyright (C) 2020 Hiroshi C. Ito
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 2 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/

#' Simulator of adaptive evolution.
#' @aliases simevol simevol-package
#' @keywords internal
"_PACKAGE"

library(deSolve);
##library(envstocker);

#' @export
getzi <- function(z,i){
    z1=z;
    for(j in 1:length(z1))z1[[j]]=(z[[j]])[i];
    return(z1);
    ##lapply(z,function(z1)(z1[i]));
}

addz <- function(z,z1){
    for(j in 1:length(z)){
        z[[j]]=c(z[[j]],z1[[j]]);
    }
    return(z);
}

max_fitness <- function(z,en){
    buf=0.0;
    for(i in 1:length(z[[1]])){
        buf1=abs(fitness(getzi(z,i),z,en));
        if(buf<buf1)buf=buf1;
    }

    return(buf);

}

#' @export
entercom <- function(){
    coms=readline("enter command: ");
    if(coms!=""){
        ##print(try(eval(parse(text=coms),envir=.GlobalEnv )));
        print(try(eval(parse(text=coms),envir=.GlobalEnv)));
        entercom();
        }
    else{
        cat("command-mode ended\n");
    }
    
}

mydir=".simevol/";
file_command=paste0(mydir,"command.R");
file_outcount=paste0(mydir,"outcount.dat");
#' @export
plot_var<-function(xid=1,yid=2){
    a$param$xid<<-xid;
    a$param$yid<<-yid;
    plot_func();
}

#' @export
pvar <- function(xid=1,yid=2){
   cat(sprintf("plot_var(%d,%d)",xid,yid),file=file_command);
}


#' @export
com <- function(void){
    cat("entercom()\n",file=file_command);
}

#' @export
halt <- function(void){
    print("simulation ended.");
    a$flag_halt<<-TRUE;
}


#' @export
hal <- function(void){
    cat("halt()\n",file=file_command);
}

#' @export
resetrange<-function(){
    a$param$xlim<<-c(NULL);
    a$param$ylim<<-c(NULL);
    plot_func();
}

#' @export
setrr<-function(){
   cat("resetrange()",file=file_command);
}

#' @export
setrange <-function(x0=NULL,x1=NULL,y0=NULL,y1=NULL){
    
    a$param$xlim<<-c(x0,x1);
    a$param$ylim<<-c(y0,y1);
        
    plot_func();
    
}

#' @export
setr <- function(x0=NULL,x1=NULL,y0=NULL,y1=NULL){
   cat(sprintf("setrange(%f,%f,%f,%f)",x0,x1,y0,y1),file=file_command);
}

#' @export
setrx <- function(x0=NULL,x1=NULL){
   cat(sprintf("setrange(%f,%f,,)",x0,x1),file=file_command);
}

#' @export
setry <- function(y0=NULL,y1=NULL){
   cat(sprintf("setrange(,,%f,%f)",y0,y1),file=file_command);
}

remove_dead1 <- function(m_next,zb,edge_die=NULL){
    die=c(NULL);
    mspe=length(zb[[1]]);
    m=as.numeric(m_next[2,]);
    m=m[2:length(m)];    
    for(i in 1:mspe){
        if(length(which(m[i]<edge_die))>0){
            die=c(die,i);
        }
    }

    
    
    if(length(die)>0){
        m=m[-die];
        for(i in 1:length(zb)){
            val=zb[[i]];
            val=val[-die];
            zb[[i]]=val;
        }
        
    }

    return(list(z=zb,en=m));
}

#' @export
invader <- function(z,en,timen=NULL,amp_invf=NULL){
    n=en[1:length(z[[1]])];
    sum_n=sum(n);
    flag_invade=0;
    x_mut=0.0;
    nsp=length(n);
    while(flag_invade==0){
    timen = timen-amp_invf*1.0/sum_n*log(runif(1)+1e-20);
    buf=0.0;
    m_spe=0;
    rv=runif(1);
    while((rv >= buf)&&(m_spe<nsp)){
    m_spe=m_spe+1;
    buf=buf+n[m_spe]/sum_n;
    }
    zpar=getzi(z,m_spe);
    mutant=mutate(zpar);
    ##fb=fitness(mutant,z,n);
    
    fb=fitness(mutant,z,en)-fitness(zpar,z,en);
    
    if(fb*amp_invf>1.0){
        cat("fit_over!! ",fb,fb*amp_invf,"\n");
        a$fit_over <<- c(a$fit_over,fb);
        }
    if( fb*amp_invf> runif(1)){
        flag_invade=1;
    }
}

return(list(z=mutant,f=fb,timen=timen));
}

output <- function(timen,nspe,z,en,file_data){
    n=en[1:length(z[[1]])];    
    cat(timen,nspe,length(z)+1,"\n",file=file_data,append=TRUE);
    for(i in 1:length(z))cat(z[[i]],"\n",file=file_data,append=TRUE);
    cat(n,"\n",file=file_data,append=TRUE);
    cat("\n",file=file_data,append=TRUE);
}


#' @export
cur.dev <- function(){
    v=as.numeric(dev.list()[length(dev.list())]);
    if(length(v)==0)v=1;
return(v);
}

#' @export
pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=100*100.0/density,outeps=FALSE,prefix="./",show=TRUE,outpng=TRUE){
    
    dens=as.character(density);
    geom=as.character(as.integer(geometry));
    dev.set(dev_id);

    tmpf=paste0(mydir,"R_pngout_temp");
    
    dev.copy2eps(file=sprintf("%s.eps",tmpf));
   
     if(outeps){
         ##system(paste("cp .R_pngout_temp.eps ",prefix,plotfile,".eps",sep=""));
         system(sprintf("cp %s.eps %s%s.eps",tmpf,prefix,plotfile));
         cat("eps output:",paste(prefix,plotfile,".eps",sep=""),"\n")
     }
    
    if(outpng){
        system(sprintf("convert -density %sx%s -geometry %s%%  -background white -alpha remove %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        system(sprintf("mv %s.png %s%s.png",tmpf,prefix,plotfile));
        cat(sprintf("png output: %s%s.png \n",prefix,plotfile));
        if(show)system(sprintf("display %s%s.png&",prefix,plotfile));
        ##system(paste("convert -density ",dens,"x",dens," -geometry ", geom,"%","  -background white -alpha remove .R_pngout_temp.eps .R_pngout_temp.png",sep=""));
       ## system(paste("mv .R_pngout_temp.png ",prefix,plotfile,".png",sep=""));
        ##cat("png output:",paste(prefix,plotfile,".png",sep=""),"\n");
        
##        if(show)system(paste("display ",prefix,plotfile,".png&",sep=""));
    }
    system(sprintf("rm %s.eps",tmpf));
    
}




    
calc_land <- function(z,en,xid=1,yid=2,xmin=-1.0,xmax=1.0,ymin=-1.0,ymax=1.0,ndiv=64){
        n=en[1:length(z[[1]])];
        xx=seq(xmin,xmax,length=ndiv);  ##とり得る形質値xのセット
        yy=seq(ymin,ymax,length=ndiv);  ##とり得る形質値yのセット
        land=matrix(ndiv*ndiv,nrow=ndiv,ncol=ndiv)*0.0;
        z1=z;
        for(i in 1:ndiv){
            for(j in 1:ndiv){
                if(length(z)>2){
                    for(k in 1:length(z))z1[[k]]=sum(z1[[k]]*n)/sum(n);
                }

                z1[[xid]]=xx[i];
                z1[[yid]]=yy[j];

            
                land[i,j]=fitness(z1,z,en);
            }
        }
        return(list(x=xx,y=yy,fit=land));
    }


pop_dynamics0 <- function(t,n,parms){
    
    dn=n*0;
    for(i in 1:length(dn)){
          dn[i]=n[i]*fitness(getzi(parms$z,i),parms$z,n);
        
            

    }
    
    list(dn);
}

#' @export
do_pop <- function(z1,en1,pop_dynamics,set_parms=NULL,times=10^seq(0,20,by=1.0)-1,edge_die=1e-6,edge_fit=1e-12){
##    .ee.append("do_pop",environment())

    for(irad in 1:length(times)){
        if(class(set_parms)=="function"){parm=set_parms(z1,en1);}
        else{parm=list(z=z1);}
        n_next=deSolve::radau(func=pop_dynamics,time=times[irad:(irad+1)],y=en1,parms=parm,hini=1e-5,hmax=1e9);
    
        res=remove_dead1(n_next,z1,edge_die=edge_die);
        z1=res$z;
        en1=res$en;
        
        if(max_fitness(z1,en1)<edge_fit)break;        
    }
    return(list(z=z1,en=en1,irad=irad));
}

#' @export
plot_init <- function(edim){

        X11(width=5,height=5);
        par(pty=a$param$win_style);  
        win0=cur.dev();
        
        npanel=length(a$z)+edim;
        ncols=npanel;
        ##nrows=as.integer(npanel/ncols)+as.integer(npanel%%ncols>0);
        nrows=1;
        om_left=1;
        om_right=0;
        om_bottom=1;
        om_top=0;
    
    ##X11(width=ncols*3.0,height=nrows*4);
        X11(width=7,height=nrows*4);
        win1=cur.dev();
        par(oma = c(om_bottom, om_left, om_top, om_right),mar=c(3,3,2,3),mgp=c(2,0.7,0)); ##bottom,left,top,right
        ##par(family=font_family) ;
        par(mfrow=c(nrows,ncols));
        return(c(win0,win1));
    }

#' @export
param0=list(
    fitness_contour=TRUE,
    time_lab="time",
    cex.lab=1.1,
    xid=1,
    yid=2,
    xlim=c(NULL),
    ylim=c(NULL),
    win_style=NULL,
    fcon=list(levels=c(0.0,0.1),col=c("red","gray"),lwd=1,lty=c(1,1)),
    resi=list(col="black",bg="green",cex=0.5,pch=21,amp=0.2,ampn=0.01),
    traj=list(col="blue",cex=0.3,pch=17),
    env=list(col="orange",lwd=1)
);

adjust_n <- function(n,amp,ampn,offset){
return (offset+amp*(log(n*ampn+1)/log(ampn+1)));
}

#' @export
plot_func <- function(){

    z=a$z;
    en=a$en;
    n=a$n;
    edim=a$edim;
    zz=a$zz;
    nn=a$nn;
    tt=a$tt;
    ee=a$ee;
    te=a$te;

    param=a$param;
    nspe=length(z[[1]]);
    timen=te[length(te)];
    dev.set(a$winid[1]);
    xid=param$xid;yid=param$yid;
    trait_names=param$trait_names;
    env_names=param$env_names;
    plot_mask=param$plot_mask;
    
    if(length(trait_names)==0)trait_names=names(z);
    if((length(env_names)==0)&&(edim>0))env_names=paste0("Env",seq(edim));
        
    if(length(z)==1){
        ndiv=128;
        ranx=1.2*max(abs(zz[[1]]));
        xx1=seq(-ranx,ranx,,ndiv);
        land=xx1*0.0;
        fit1=z[[1]]*0.0;
        x1=z[[1]];
        z1=getzi(z,1);
        for(i in 1:length(xx1)){
            z1[[1]]=xx1[i];
            land[i]=fitness(z1,z,en);
        }
        for(i in 1:length(fit1)){
            z1[[1]]=x1[i];
            fit1[i]=fitness(z1,z,en);
        }
        plot(xx1,land,col=param$fcon$col[1],lty=param$fcon$lty[1],type="l",xlab=trait_names[1],ylab="Fitness",ylim=c(-max(land)*0.2,max(land)*1.0),cex.lab=param$cex.lab);
        lines(xx1,land*0,col=param$fcon$col[2],lty=param$fcon$lty[2]);
        
        points(x1,fit1,col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
        
    }
    else{
        
        ndiv=64;


        if((yid>0)&&(xid>0)){
        rany=1.2*max(abs(zz[[yid]]));        
        ranx=1.2*max(abs(zz[[xid]]));        
        if(param$fitness_contour){
            land=calc_land(z,en,xid,yid,xmin=-ranx,xmax=ranx,ymin=-rany,ymax=rany);
        }
        
        plot(zz[[xid]],zz[[yid]],col=param$traj$col,pch=param$traj$pch,cex=param$traj$cex,xlab=trait_names[xid],ylab=trait_names[yid],cex.lab=param$cex.lab,xlim=param$xlim,ylim=param$ylim);
        if(param$fitness_contour)contour(land$x,land$y,land$fit,add=1,levels=param$fcon$levels,col=param$fcon$col,lty=param$fcon$lty);
        
        points(z[[xid]],z[[yid]],col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
        }
        else{
            if(xid>0){
            plot(zz[[xid]],tt,col=param$traj$col,pch=param$traj$pch,cex=param$traj$cex,xlab=trait_names[xid],ylab="Time",cex.lab=param$cex.lab,xlim=param$xlim,ylim=param$ylim);
            points(z[[xid]],rep(timen,nspe),col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
            }
            if(yid>0){
                plot(tt,zz[[yid]],col=param$traj$col,pch=param$traj$pch,cex=param$traj$cex,ylab=trait_names[yid],xlab="Time",cex.lab=param$cex.lab,xlim=param$xlim,ylim=param$ylim);
                points(rep(timen,nspe),z[[yid]],col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
            }


            
        }


        
    }

    
    dev.set(a$winid[2]);
    for(i in 1:length(z)){
        if(class(plot_mask[[i]])=="function"){
            lis=which(plot_mask[[i]](z)==1);
            liss=which(plot_mask[[i]](zz)==1);
            plot((zz[[i]])[liss],tt[liss],col=param$traj$col,pch=param$traj$pch,cex=param$traj$cex,xlab=trait_names[i],ylab="Time",cex.lab=param$cex.lab);
            points((z[[i]])[lis],rep(timen,nspe)[lis],col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
        }
        else{            
            plot(zz[[i]],tt,col=param$traj$col,pch=param$traj$pch,cex=param$traj$cex,xlab=trait_names[i],ylab="Time",cex.lab=param$cex.lab);
            points(z[[i]],rep(timen,nspe),col=param$resi$col,bg=param$resi$bg,pch=param$resi$pch,cex=adjust_n(n,param$resi$amp,param$resi$ampn,param$resi$cex));
      

        }
        }
    if(edim>0){
        for(i in 1:edim){
        plot(ee[,i],te,col=param$env$col,lwd=param$env$lwd,type="l",xlab=env_names[i],ylab="Time",cex.lab=param$cex.lab);  
    }
    }
}

#' @export
comwin <- function(make_win=0){
    if(0){
    defs=paste('file_command="command.R";',
               'setr <- function(x0,x1,y0,y1){',
               'cat(sprintf("setrange(%f,%f,%f,%f)",x0,x1,y0,y1),file=file_command);',
               '}',
               ' ',
               'setrx <- function(x0,x1){',
               'cat(sprintf("setrange(%f,%f,,)",x0,x1),file=file_command);',
               '}',
               ' ',
               'setry <- function(y0,y1){',
               'cat(sprintf("setrange(,,%f,%f)",y0,y1),file=file_command);',
               '}',
               ' ',
               'setrr <- function(){',
               'cat("resetrange()",file=file_command);',
               '}',
               ' ',
               'com <- function(void){',
               'cat("entercom()",file=file_command);',
               '}',
               ' ',
               'hal <- function(void){',
               'cat("halt()",file=file_command);',
               '}',
               ' ',
               sep="\n");

    cat(defs,file=".Rprofile");
    ##Sys.sleep(0.1);
    }
    cat("library(simevol);\n",file=".Rprofile");
    if(make_win!=0)system("gnome-terminal --geometry=40x10 -- R");
    
}

##' function for evolutionary simulation.
##' @title function for evolutionary simulation
##' @param z : traits
##' @author Hiroshi C. Ito
#' @export
simevol <- function(z,en,fitness=NULL,mutate=NULL,pop_dynamics=NULL,set_parms=NULL,
                  edge_die=1e-6,
                  edge_fit=1e-12,
                  tmax=100000000,
                  times=10^seq(0,20,by=1.0)-1,
                  out_interval=10,
                  show_interval=10,
                  amp_invf=0.1,
                  level_invf=0.01,
                  file_data="test.dat",
                  fitness_contour=TRUE,
                  yorick_plot=NULL,
                  plot_mask=NULL,
                  reset_win=TRUE,
                  trait_names=NULL,
                  env_names=NULL,
                  param=param0){

    ##    .ee.append("adsim",environment())

   ##if(length(trait_names)==0)trait_names=names(z);

 
    
    if(class(pop_dynamics)!="function")pop_dynamics=pop_dynamics0;    
    timen=0.0;

    outcount=1;

    nspe=length(z[[1]]);
    n=en[1:nspe];

    edim=length(en)-nspe;

    param=c(param,list(plot_mask=plot_mask,trait_names=trait_names,env_names=env_names));
    
    a<<-list(
        z=z,
        en=en,
        n=n,
        edim=edim,
        zz=z,
        nn=c(n),
        tt=c(rep(0.0,nspe)),
        ee=c(NULL),
        te=c(0.0),
        winid=c(2,3),
        invf=c(0.1),
        fit_over=c(NULL),
        param=param,
        flag_halt=FALSE
        );
    
    if(a$edim>0)a$ee<<-c(a$ee,en[(nspe+1):length(en)]);
    
    if(file.exists(file_data))file.remove(file_data);

    output(timen,length(n),z,n,file_data);
    cat(outcount,file=file_outcount);
    outcount=outcount+1;
    
    if(class(yorick_plot)=="function"){
        yorick_plot();
        system("xterm -fn 7x14 -bg navy -fg white -e 'rlwrap -c ./yorick_idl_follow.sh'&");
    }
    
    if(reset_win || (length(dev.list())<2)){
        graphics.off();
        a$winid<<-plot_init(edim);
    }
    res=do_pop(z,en,pop_dynamics,set_parms=set_parms,times=times,edge_die=edge_die,edge_fit=edge_fit);
    z=res$z;
    nspe=length(z[[1]]);
    en=res$en;
    n=en[1:nspe];


    if(!file.exists(mydir))system(sprintf("mkdir %s",mydir));
    
    if(file.exists(file_command))file.remove(file_command);

   
    for(t in 1:tmax){
        if(a$flag_halt==T)break;
        if(file.exists(file_command)){
            source(file_command);
            file.remove(file_command);
        }
        
    inv=invader(z,en,timen,amp_invf);
    timen=inv$timen;

        
    amp_invf=max(level_invf,level_invf/mean(a$invf[max((t-100),1):length(a$invf)]));
    ##aa=c(aa,amp_invf);

        nspe=length(z[[1]]);
        z=addz(z,inv$z);
        a$invf<<-c(a$invf,inv$f);    
        n=c(n,edge_die*10);
        
        if(edim>0){en=c(n,en[(nspe+1):length(en)]);}
        else{en=n;}
       
        nspe=length(n);
        
        res=do_pop(z,en,pop_dynamics,set_parms=set_parms,times=times,edge_die=edge_die,edge_fit=edge_fit);
       
        z=res$z;
        nspe=length(z[[1]]);
        en=res$en;
        n=en[1:nspe];
        
        a$z<<-z;
        a$en<<-en;
        a$n<<-n;

        if(t%%out_interval==0){
            output(timen,nspe,z,n,file_data);
            cat(outcount,file=file_outcount);
            outcount=outcount+1;
    a$zz<<-addz(a$zz,z);
    if(a$edim>0)a$ee<<-rbind(a$ee,en[(nspe+1):length(en)]);
            a$nn<<-c(a$nn,n);
            a$tt<<-c(a$tt,rep(timen,nspe));
            a$te<<-c(a$te,timen);
            
}
    
    if(t%%show_interval==0){
    cat("time:",timen, " residents:",nspe, " invasion:", t,"amp:",amp_invf,"fit_over:",length(a$fit_over),"irad:",res$irad,"\n");

param0$fitness_contour=fitness_contour;
    plot_func();    
    
    }
}

}


