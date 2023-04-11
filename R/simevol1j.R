# file simevol/R/simevol1f.R
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

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/      Functions for simulation      _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

#' @export
geti <- function(phe,i){
    phe1=phe;
    for(j in 1:length(phe1))phe1[[j]]=(phe[[j]])[i];
    return(phe1);
    ##lapply(z,function(z1)(z1[i]));
}

#' @export
add_phenotype <- function(phe,phe1){
    for(j in 1:length(phe)){
        phe[[j]]=c(phe[[j]],phe1[[j]]);
    }
    return(phe);
}

#' @export
fitness_all <- function(phe,en){
    buf=phe[[1]]*0.0;
    for(i in 1:length(phe[[1]]))buf[i]=.simevol_func$fitness(geti(phe,i),phe,en);
    return(buf);
}


max_fitness <- function(phe,en){
    buf=0.0;
    for(i in 1:length(phe[[1]])){
        buf1=abs(.simevol_func$fitness(geti(phe,i),phe,en));
        if(buf<buf1)buf=buf1;
    }
    return(buf);
}


#' @export
pop_dynamics0 <- function(t,n,parms){    
    dn=n*0;
    for(i in 1:length(dn)){
          dn[i]=n[i]*.simevol_func$fitness(geti(parms$phe,i),parms$phe,n);
    }
    list(dn);
}

#' @export
mutate0 <- function(phe){
    mut=phe;
    for(i in 1:(length(phe)-2))mut[[i]]=phe[[i]]+rnorm(1,mean=0.0,sd=a$sparam$m_sd);
    return(mut);
}


#' @export
invader <- function(phe=a$phe,en=a$en,mutate=.simevol_func$mutate,timen=0.0,amp_invf=a$sparam$amp_invf){
    n=en[1:length(phe[[1]])];
    sum_n=sum(n);
    flag_invade=0;
    x_mut=0.0;
    nsp=length(n);
    while(flag_invade==0){
        timen = timen-(1.0/a$sparam$m_rate)*amp_invf*(1.0/sum_n)*log(runif(1)+1e-20);
        buf=0.0;
        m_spe=0;
        rv=runif(1);
        while((rv >= buf)&&(m_spe<nsp)){
            m_spe=m_spe+1;
            buf=buf+n[m_spe]/sum_n;
        }
        phe_par=geti(phe,m_spe);
        mutant=mutate(phe_par);
        
        ##fb=fitness(mutant,phe,en)-fitness(phe_par,phe,en);
        fb=fitness(mutant,phe,en);
        
        if(fb*amp_invf>1.0){
            cat("Invasion fitness larger than 1 !! ",fb,fb*amp_invf,"\n");
            a$sparam$fit_over <<- c(a$sparam$fit_over,fb);
        }
        if( fb*amp_invf> runif(1)){
            flag_invade=1;
        }
    }
    mutant=c(mutant,list(t=timen,pid=a$ninv+1));
    return(list(phe=mutant,n=a$sparam$n_mutant_init,f=fb,timen=timen,pid_par=phe_par$pid));
}



#' @export
add_invader <- function(phe,en,inv){    
    edim=length(en)-length(phe[[1]]);
    nspe=length(phe[[1]]);
    ebuf=c(NULL);
    if(edim>0)ebuf=en[(nspe+1):length(en)];
    n=en[1:nspe];

    phe=add_phenotype(phe,inv$phe);        
    n=c(n,inv$n);
    en=c(n,ebuf);
    nspe=length(n);

    return(list(phe=phe,en=en,n=n,nspe=nspe,inv=inv));
}



#' @export
remove_extinct <- function(en,phe,edge_extinct=NULL){
    die=c(NULL);
    nspe=length(phe[[1]]);
    
    for(i in 1:nspe){
        if(length(which(en[i]<edge_extinct))>0){
            die=c(die,i);
        }
    }

    
    if(length(die)>0){
                
        en=en[-die];
        for(i in 1:length(phe)){
            val=phe[[i]];
            val=val[-die];
            phe[[i]]=val;
        }
        
    }

    return(list(phe=phe,en=en,die=die));
}


#' @export
rootfunc <- function(t,n,parms){
    if(parms$edim==0){
        return(n-0.1*parms$edge_extinct);
    }
    else{
        n=n-0.1*parms$edge_extinct;
        n[(parms$nspe+1):length(n)]=1.0;
        return(n);
    }
}

#' @export
eventfunc <- function(t,n,parms){
    lis=which(n<parms$edge_extinct);
     if(length(lis)>0)n[lis]=n[lis]*0.0;
     
    return(n);
}


#' @export
simpop_check <- function(phe=a$phe,en=a$en,inv=invader(),pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,check=TRUE,nrad=7,drad=1.0,divrad=64,xlim=c(NULL),out=FALSE,reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
    res=simpop_invade(phe=phe,en=en,inv=inv,pop_dynamics=pop_dynamics,set_parms=set_parms,edge_extinct=edge_extinct,edge_fit=edge_fit,check=check,nrad=nrad,drad=drad,divrad=divrad,xlim=xlim,reset_win=reset_win,logt=logt,param_desolve=param_desolve);
    if(out==TRUE)return(res);
}
    
#' @export
simpop_invade <- function(phe=a$phe,en=a$en,inv=a$inv,pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,check=FALSE,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad,xlim=c(NULL),reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
##.ee.append("simpop_invade",environment())

    state=add_invader(phe,en,inv);       
    state_new=simpop(state$phe,state$en,pop_dynamics=pop_dynamics,set_parms=set_parms,edge_extinct=edge_extinct,edge_fit=edge_fit,check=check,nrad=nrad,drad=drad,divrad=divrad,xlim=xlim,reset_win=reset_win,logt=logt,param_desolve=param_desolve);
    return(state_new);
    
}


#' @export
simpop_set_parms <- function(phe1,en1,set_parms=NULL){
            if(class(set_parms)=="function"){
                parm=set_parms(phe1,en1);
            }else{
                parm=list(phe=phe1);               
            }
return(parm);
}

#' @export
get_last_en <- function(en_next){
    lent=length(en_next[,1]);
    en=as.numeric(en_next[lent,]);
    return(en[2:length(en)]);    
}


#' @export
simpop <- function(phe1=a$phe,en1=a$en,pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extince,edge_fit=a$sparam$edge_fit,check=FALSE,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad,xlim=c(NULL),reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
    hini=param_desolve$hini;
    hmax=param_desolve$hmax;
    rtol=param_desolve$rtol;
    atol=param_desolve$atol;
    method=param_desolve$method;
    
##.ee.append("simpop",environment())
    nspe0=length(phe1[[1]]);
    nspe1=nspe0;
    edim1=length(en1)-nspe1;
    if(check==FALSE){
        for(irad in 1:nrad){            
            parm=c(simpop_set_parms(phe1,en1,set_parms),list(edge_extinct=edge_extinct,nspe=nspe1,edim=edim1));
             mytime=10^seq(drad*irad,drad*(irad+1),,divrad);
            n_next=deSolve::ode(y=en1,times=mytime,func=pop_dynamics,parms=parm,method=method,hini=hini,hmax=hmax,rootfun=rootfunc,rtol=rtol,atol=atol);
            
            res=remove_extinct(get_last_en(n_next),phe1,edge_extinct=edge_extinct);
            phe1=res$phe;
            en1=res$en;
            nspe1=length(phe1[[1]]);
            edim1=length(en1)-nspe1;
            max_fit=max(abs(fitness_all(phe1,en1)));
            max_t=max(n_next[,1]);
            ##if((max_t<times[irad+1])&&(length(phe1[[1]])==nspe0)){
            ##    cat("population dynamis failed!! time:", max_t,"max_fit:",max_fit,"\n");
            ##    a$sparam$flag_halt<<-T;
            ##}

            if(length(which(res$die==nspe0))>0){
                cat("inveder extinct!!\n");
                ##if(a$sparam$error_stop==TRUE);                
                ##a$sparam$flag_halt<<-T;
                }
                        

            if(max_fit<edge_fit)break;
            ##hini=min(10^(drad*irad),1e-3/max_fit);
            hini=min(10^(drad*irad-1),1e-4/max_fit); ## this could be improved
            ##cat("hini:",hini,"\n")
        }
        
        return(list(phe=phe1,en=en1,irad=irad));
        
    }
    
    if(check==TRUE){
        phe0=phe1;
        en0=en1;
        
        parm=c(simpop_set_parms(phe1,en1,set_parms),list(edge_extinct=edge_extinct,nspe=nspe1,edim=edim1));        
        for(irad in 1:nrad){    
            mytime=10^seq(drad*irad,drad*(irad+1),,divrad);
            n_next=deSolve::ode(y=en1,times=mytime,func=pop_dynamics,parms=parm,method=method,hini=hini,hmax=hmax,rootfun=rootfunc,atol=atol,rtol=rtol);
            en1=(n_next[nrow(n_next),])[2:ncol(n_next)];
            if(irad==1)n_nex=n_next;
            if(irad>1)n_nex=rbind(n_nex,n_next);
            max_fit=max(abs(fitness_all(phe1,en1)));
            if(max_fit<edge_fit)break;
            hini=min(10^(drad*irad-1),1e-4/max_fit);
        }
        n_next=n_nex;
      
        res=remove_extinct(get_last_en(n_next),phe1,edge_extinct=edge_extinct);
        phe1=res$phe;
        en1=res$en;
        if(length(which(res$die==nspe0))>0){
                cat("inveder extinct!!\n");
        }
        plot_simpop(n_next,phe0,edge_extinct,reset_win=reset_win,logt=logt,xlim=xlim);
        print("fitness");
        print(fitness_all(phe1,en1));
        
        return(list(phe=phe1,en=en1,n_next=n_next));
        
        
    }
}

#' @export
plot_simpop <- function(n_next,phe0,edge_extinct,reset_win=FALSE,logt=FALSE,xlim=NULL){
    nspe=length(phe0[[1]]);
    edim=ncol(n_next)-nspe-1;
    
            nam=colnames(n_next);
            name_n=nam[2:(nspe+1)];
            colnames(n_next)[2:(nspe+1)]=paste0("n",name_n);
            
            en_last=n_next[nrow(n_next),];
            n_last=en_last[2:(nspe+1)];
            t_last=en_last[1];
            if(edim>0){
                name_e=nam[(nspe+2):length(nam)];
                colnames(n_next)[(nspe+2):length(nam)]=paste0("e",seq(1,length(name_e)));
                e_last=en_last[(nspe+2):length(en_last)];
            }
            lis=which(n_last<edge_extinct);
            color=rep("blue",ncol(n_next)-1);
            color[nspe]="green";
            if(length(lis)>0)color[lis]="red";
            
            npanel=nspe+edim;
            ncols=3;
            nrows=as.integer(npanel/ncols)+as.integer(npanel%%ncols>0);
            
        if((reset_win==FALSE)&&(a$winid[3]>0)){
            dev.set(a$winid[3]);
            plot.new();
        }
        else{
            X11(width=ncols*3.0,height=nrows*1.5);
            a$winid[3]<<-cur.dev();
        }
        om_left=1;
        om_right=0;
        om_bottom=1;
        om_top=0;

        par(oma = c(om_bottom, om_left, om_top, om_right),mar=c(3,3,2,3),mgp=c(2,0.7,0)); ##bottom,left,top,right
       
        par(mfrow=c(nrows,ncols));
        name=colnames(n_next);
        for(i in 1:(ncol(n_next)-1)){
            if(i<=nspe){
                ylim=c(0,max(n_next[,(i+1)])*1.2);
                ##ylim=c(NULL);
                tit="";
               
                for(j in 1:(length(phe0)-2))tit=paste(tit,sprintf("%s:%f",names(phe0)[j],phe0[[j]][i]));
                if(logt){plot(n_next[,1],n_next[,(i+1)],log="x",col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,xlim=xlim,main=tit);}
                else{plot(n_next[,1],n_next[,(i+1)],col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,xlim=xlim,main=tit);}
            }
            else{
                if(logt){plot(n_next[,1],n_next[,(i+1)],log="x",col="orange",type="l",xlab=name[i+1],ylab="value",xlim=xlim);}
                else{plot(n_next[,1],n_next[,(i+1)],col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,main=tit,xlim=xlim);}
            }
        }
        print(phe0);
        print(t_last);
        print(n_last);
        if(edim>0)print(e_last);
}


##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/        Plotting functions          _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

#' @export
calc_fitness_land <- function(phe,en,xid=1,yid=2,xmin=-1.0,xmax=1.0,ymin=-1.0,ymax=1.0,ndiv=64){
    n=en[1:length(phe[[1]])];
   
    xx=seq(xmin,xmax,,ndiv);
    yy=seq(ymin,ymax,,ndiv); 
    
    land=matrix(ndiv*ndiv,nrow=ndiv,ncol=ndiv)*0.0;

    if(length(a$pparam$fitness_contour_phe)>0){
        phe1=a$pparam$fitness_contour_phe;
    }
    else{
        phe1=phe;
        pdim=length(phe)-2;
        if(pdim>1){
            for(k in 1:pdim)phe1[[k]]=sum(phe1[[k]]*n)/sum(n);
        }
    }
    for(i in 1:ndiv){
        for(j in 1:ndiv){
            phe1[[xid]]=xx[i];
            phe1[[yid]]=yy[j];
            land[i,j]=.simevol_func$fitness(phe1,phe,en);
        }
    }
    return(list(x=xx,y=yy,fit=land));
}

#' @export
plot_fitness_contour <- function(land,p){
    contour(land$x,land$y,land$fit,add=1,levels=p$fcon$levels,col=p$fcon$col,lty=p$fcon$lty);
}


#' @export
plot_init <- function(reset_win=TRUE){
    if(reset_win){
    X11(width=5,height=5,title=paste("runid:",a$runid,"  ",a$runname,"    dev:",(cur.dev()+1)));
    a$winid[1] <<-cur.dev();
    }
    else{
        dev.set(a$winid[1]);
    }
    par(pty=a$pparam$win_style);
        
    npanel=length(a$phe)-2+a$edim;
    ncols=npanel;
    ##nrows=as.integer(npanel/ncols)+as.integer(npanel%%ncols>0);
    nrows=1;
    om_left=1;
    om_right=0;
    om_bottom=1;
    om_top=0;

   
    if(a$show_subwin==TRUE){
        if(reset_win){
        ##X11(width=ncols*3.0,height=nrows*4);
            X11(width=7,height=nrows*4,title=paste("runid:",a$runid,"  ",a$runname,"    dev:",(cur.dev()+1)));

            a$winid[2]<<-cur.dev();
        }
        else{
            dev.set(a$winid[2]);        
        }
        
        par(oma = c(om_bottom, om_left, om_top, om_right),mar=c(3,3,2,3),mgp=c(2,0.7,0)); ##bottom,left,top,right
        ##par(family=font_family) ;
        par(mfrow=c(nrows,ncols));
    }
        ##return(c(win0,win1));
    }

#' @export
pparam0=list(
    fitness_contour=TRUE,
    fitness_contour_phe=NULL,
    time_lab="Time",
    cex.lab=1.1,
    xid=1,
    yid=2,
    xlim=c(NULL),
    ylim=c(NULL),
    win_style=NULL,
    fcon=list(levels=c(0.0,0.1),col=c("red","gray"),lwd=1,lty=c(1,1)),
    resi=list(col="black",bg="green",cex=0.4,pch=21,amp=1.0,ampn=10,bg2="red",bgid=-1),
    traj=list(col="blue",cex=1.0,pch=17,every=10,lwd=0.5),
    env=list(col="orange",lwd=1),
    plot_mask=NULL,
    trait_names=NULL,
    nv_names=NULL,
    fitness_contour=TRUE,
    fitness_contour_phe=NULL,
    pal=(rev(rainbow(100,end=0.7))),
    palfunc=NULL
);

#' @export
adjust_n <- function(n,amp,ampn,offset){
    return (offset+amp*(log(n*ampn+1)/log(ampn+1)));
}

#' @export
plot_lim <- function(xid,yid,p){
    if((xid>0) && (yid>0)){
        xp=c(min(a$traj$phe[[xid]]),max(a$traj$phe[[xid]]));
        yp=c(min(a$traj$phe[[yid]]),max(a$traj$phe[[yid]]));
        xlab=p$trait_names[xid];
        ylab=p$trait_names[yid];
    }
    else{
        if(xid>0){
            xp=c(min(a$traj$phe[[xid]]),max(a$traj$phe[[xid]]));
            yp=c(min(a$traj$t),max(a$traj$t));
            xlab=p$trait_names[xid];
            ylab=p$time_lab;
            
        }
        if(yid>0){
            xp=c(min(a$traj$t),max(a$traj$t));
            yp=c(min(a$traj$phe[[yid]]),max(a$traj$phe[[yid]]));
            xlab=p$time_lab;
            ylab=p$trait_names[yid];
            
        }
        
    }
    

        plot(xp,yp,type="n",,cex.lab=p$cex.lab,xlim=p$xlim,ylim=p$ylim,xlab=xlab,ylab=ylab);
}

#' @export
plot_lim_sub <- function(xid,p){
    xp=c(min(a$traj$phe[[xid]]),max(a$traj$phe[[xid]]));
    yp=c(min(a$traj$t),max(a$traj$t));
    xlab=p$trait_names[xid];
    ylab=p$time_lab;
    
    plot(xp,yp,type="n",,cex.lab=p$cex.lab,xlim=p$sub_xlim[[xid]],ylim=p$sub_ylim[[xid]],xlab=xlab,ylab=ylab);
}


adj_tree <- function(tree){
    for(i in 1:length(tree)){
        b=tree[[i]];
        lis=which(a$phe$pid==b$pid[length(b$pid)]);
        if(length(lis)>0){
            b=geti(b,length(b$pid));
            b$t=a$timen;
            tree[[i]]=add_phenotype(tree[[i]],b);
        }
    }
    return(tree);
}

        
pline <- function(b,xid,yid,...){
        lines(b[[xid]],b[[yid]],...);
}

#' @export
stree <- function(ename){
    lapply(a$tree,eval(parse(text=sprintf('function(x)x$%s',ename))));
}



#' @export
plot_traj_line <-function(tree,xid,yid,p,mask=NULL){
    if(xid==0)xid=a$pdim+1;
    if(yid==0)yid=a$pdim+1;
        lapply(tree,pline,xid,yid,col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex,lwd=p$traj$lwd);    
}


#' @export
plot_traj <-function(xid,yid,p,mask=NULL){
    if(length(mask)==0){
        if((xid>0)&&(yid>0)){
            points(a$traj$phe[[xid]],a$traj$phe[[yid]],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
        else{        
            if(xid>0)points(a$traj$phe[[xid]],a$traj$t,col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
            if(yid>0)points(a$traj$t,a$traj$phe[[yid]],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
    }
    else{
        if((xid>0)&&(yid>0)){
            points((a$traj$phe[[xid]])[mask],(a$traj$phe[[yid]])[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
        else{        
            if(xid>0)points((a$traj$phe[[xid]])[mask],(a$traj$t)[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
            if(yid>0)points((a$traj$t)[mask],(a$traj$phe[[yid]])[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
    }
}

#' @export
plot_phe <-function(xid,yid,p){
    timen=a$timen;
    nspe=length(a$phe[[1]]);
    if((xid>0)&&(yid>0)){
        points(a$phe[[xid]],a$phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    else{
        if(xid>0)points(a$phe[[xid]],rep(timen,nspe),col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
        if(yid>0)points(rep(timen,nspe),a$phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    
}

#' @export
plot_1dim <- function(p){
    phe=a$phe;
    en=a$en;
    n=a$n;
    edim=a$edim;
    nspe=length(phe[[1]]);
    
    ndiv=128;
    ranx=1.2*max(abs(a$traj$phe[[1]]));
    xx1=seq(-ranx,ranx,,ndiv);
    land=xx1*0.0;
    fit1=phe[[1]]*0.0;
    x1=phe[[1]];
    phe1=geti(phe,1);
        for(i in 1:length(xx1)){
            phe1[[1]]=xx1[i];
            land[i]=fitness(phe1,phe,en);
        }
        for(i in 1:length(fit1)){
            phe1[[1]]=x1[i];
            fit1[i]=fitness(phe1,phe,en);
        }
        plot(xx1,land,col=p$fcon$col[1],lty=p$fcon$lty[1],type="l",xlab=p$trait_names[1],ylab="Fitness",ylim=c(-max(land)*0.2,max(land)*1.0),cex.lab=p$cex.lab);
        lines(xx1,land*0,col=p$fcon$col[2],lty=p$fcon$lty[2]);
        
        points(x1,fit1,col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
}


#' @export
plot_func0 <- function(traj_line=TRUE){
    phe=a$phe;
    en=a$en;
    n=a$n;
    edim=a$edim;
    traj=a$traj;
    p=a$pparam;
    nspe=length(phe[[1]]);
    timen=a$timen;
    
    dev.set(a$winid[1]);
    
    xid=p$xid;yid=p$yid;
    trait_names=p$trait_names;
    env_names=p$env_names;
    plot_mask=p$plot_mask;

    dev.hold();
    if(class(.simevol_func$palfunc)=="function"){
            npal=length(p$pal);
            bgval=.simevol_func$palfunc(phe,en);
        
            ##cmax=max(bgval)+1e-10;
            ##cmin=min(bgval)-1e-10;
            ##cid=as.integer((npal-1)*(bgval-cmin)/(cmax-cmin))+1;
            cid=as.integer((npal-1)*bgval)+1;
            p$resi$bg=p$pal[cid];
    }else{
        if(p$resi$bgid>-1){
            bgid=p$resi$bgid;
            npal=length(p$pal);
            if(bgid==0){
                bgval=n;
            }else{
                bgval=phe[[bgid]];
            }
            cmax=max(bgval)+1e-10;
            cmin=min(bgval)-1e-10;
            cid=as.integer((npal-1)*(bgval-cmin)/(cmax-cmin))+1;
            p$resi$bg=p$pal[cid];
        }
    }
        
    
    if(length(trait_names)==0)trait_names=names(phe);
    if((length(env_names)==0)&&(edim>0))env_names=paste0("Env",seq(edim));
        
    if((length(phe)-2)==1){
        plot_1dim(p);
    }else{        
       plot_lim(xid,yid,p);
       if((xid>0)&&(yid>0)&&p$fitness_contour){
           ranx=max(1.2*max(abs(traj$phe[[xid]])),0.001);        
           rany=max(1.2*max(abs(traj$phe[[yid]])),0.001);        
           
           land=calc_fitness_land(phe,en,xid,yid,xmin=-ranx,xmax=ranx,ymin=-rany,ymax=rany);
       }            
       if(traj_line==TRUE){
           tree1=adj_tree(a$tree);
           plot_traj_line(tree1,xid,yid,p);
       }
       else{
           plot_traj(xid,yid,p);
       }
        if((xid>0)&&(yid>0)&&p$fitness_contour)plot_fitness_contour(land,p);
        plot_phe(xid,yid,p);
    }

        
    dev.flush();
    
    if(a$show_subwin==TRUE){
        dev.set(a$winid[2]);
        dev.hold();
    for(i in 1:(length(phe)-2)){
        ##plot_lim(i,0,p);
        plot_lim_sub(i,p);
        if(traj_line==TRUE){
            if((length(phe)-2)==1)tree1=adj_tree(a$tree);
                       plot_traj_line(tree1,i,0,p);
       }
       else{
                       plot_traj(i,0,p);
       }


            plot_phe(i,0,p);
        
    }
    if(edim>0){
        for(i in 1:edim){
            plot(traj$e[,i],traj$te,col=p$env$col,lwd=p$env$lwd,type="l",xlab=env_names[i],ylab="Time",cex.lab=p$cex.lab);  
        }
    }
        dev.flush();
    }

    
}

#' @export
ptraj <- function(li=1){
    if(li==1)plot_func0(traj_line=TRUE);
    if(li==0)plot_func0(traj_line=FALSE);

}

#' @export
output0 <- function(timen,z,n){
    options(scipen=100); 
    nspe=length(z[[1]]);
    cat(timen,nspe,"\n",file=a$file_data,append=TRUE);
    cat(a$phe$pid+1,"\n",file=a$file_data,append=TRUE);
    for(i in 1:nspe)cat(z$x[i],z$y[i],n[i],"\n",file=a$file_data,append=TRUE);
    cat("\n",file=a$file_data,append=TRUE);
    options(scipen=0);
}


#' @export
output_tree0 <- function(fname){
    options(scipen=100);
    q=unlist(lapply(a$tree,function(z)(length(z$pid))));
    cat(length(a$tree),max(q),"\n\n",file=fname,append=FALSE);
    for(i in 1:length(a$tree)){
        time_ed=-1.0;
        blen=length(a$tree[[i]]$pid);
        if(a$tree[[i]]$pid[blen]==-1)time_ed=a$tree[[i]]$t[blen];
        time_st=a$tree[[i]]$t[1];
        
        cat(i,blen,time_st,time_ed,"\n",file=fname,append=TRUE);
        write((a$tree[[i]]$pid+1),ncolumns=50,file=fname,append=TRUE);
        cat("\n",file=fname,append=TRUE);
    }

    cat("\n\n",file=fname,append=TRUE);

    cat(a$pdim+2,length(a$tree_phe$pid),"\n",file=fname,append=TRUE);

    for(i in 1:a$pdim){
        buf=unlist(a$tree_phe[[i]]);
        write(buf,file=fname,append=TRUE,ncolumns=50);
       cat("\n",file=fname,append=TRUE);    
    }

    
    write(a$tree_phe$t,file=fname,append=TRUE,ncolumns=50);
    cat("\n",file=fname,append=TRUE);    
    write((a$tree_phe$pid+1),file=fname,append=TRUE,ncolumns=50);
    cat("\n",file=fname,append=TRUE);
    options(scipen=0);
}



##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/     Functions for simulation controll      _/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

#' @export
comwin <- function(make_win=0,mydir=a$sparam$mydir){
    ##cat("library(simevol);\n",file=".Rprofile");
    ##cat(sprintf('source("simevol1g.R");file_command<<-"%scommand";\n',mydir),file=".Rprofile");
    cat(sprintf('library(simevol);file_command<<-"%scommand";\n',mydir),file=".Rprofile");
    if(make_win!=0)system("gnome-terminal --geometry=40x10 -- R");
    
}


##' @title send commands to simulators
#' @export
scom <- function(scom,runid=1){
    for(i in 1:length(runid))cat(sprintf("try(eval(parse(text=\'%s\')))\n",scom),file=paste0(file_command,runid[i],".R"));
}


#' @export
plot_var<-function(xid=1,yid=2){
    a$pparam$xid<<-xid;
    a$pparam$yid<<-yid;
    .simevol_func$plot_func();
}


##' @title change valuables for plotting
#' @export
pvar <- function(xid=1,yid=2,runid=1){
    for(i in 1:length(runid))cat(sprintf("plot_var(%d,%d)\n",xid,yid),file=paste0(file_command,runid[i],".R"));
}


#' @export
entercom <- function(){
    ##cat("enter command:");
    ##coms=scan("stdin",character(),n=1);
    coms=readline("enter command: ");
    if(coms!=""){
        print(try(eval(parse(text=coms),envir=.GlobalEnv)));
        entercom();
        }
    else{
        cat("command-mode ended\n");
    }
    
}


##' @title command prompt
#' @export
com <- function(runid=1){
    cat("entercom()\n",file=paste0(file_command,runid,".R"));
}


#' @export
halt <- function(void){
    print("simulation ended.");
    a$sparam$flag_halt<<-TRUE;
}


##' @title halt simulation
#' @export
hal <- function(runid=1){
    for(i in 1:length(runid))cat("halt()\n",file=paste0(file_command,runid[i],".R"));
}

#' @export
resetrange<-function(){
    a$pparam$xlim<<-c(NULL);
    a$pparam$ylim<<-c(NULL);
    .simevol_func$plot_func();
}

#' @export
setrr<-function(runid=1){
   for(i in 1:length(runid))cat("resetrange()",file=paste0(file_command,runid[i],".R"));
}



#' @export
setrange <-function(x0=NULL,x1=NULL,y0=NULL,y1=NULL){
    a$pparam$xlim<<-c(x0,x1);
    a$pparam$ylim<<-c(y0,y1);
    .simevol_func$plot_func();
}

#' @export
setrangex <-function(x0=NULL,x1=NULL){
    a$pparam$xlim<<-c(x0,x1);
    .simevol_func$plot_func();
}

#' @export
setrangey <-function(y0=NULL,y1=NULL){
    a$pparam$ylim<<-c(y0,y1);
    .simevol_func$plot_func();
}

#' @export
setr <- function(x0=NULL,x1=NULL,y0=NULL,y1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrange(%f,%f,%f,%f)",x0,x1,y0,y1),file=paste0(file_command,runid[i],".R"));
}

#' @export
setrx <- function(x0=NULL,x1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrangex(%f,%f)",x0,x1),file=paste0(file_command,runid[i],".R"));
}

#' @export
setry <- function(y0=NULL,y1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrangey(%f,%f)",y0,y1),file=paste0(file_command,runid[i],".R"));
}


#' @export
resetrange_sub<-function(id){
    a$pparam$sub_xlim[id]<<- list(NULL);
    a$pparam$sub_ylim[id]<<- list(NULL);
    .simevol_func$plot_func();
}

#' @export
setrange_sub <-function(x0=NULL,x1=NULL,id=1){
    a$pparam$sub_xlim[[id]]<<-c(x0,x1);
    .simevol_func$plot_func();
}


#' @export
cur.dev <- function(){
    v=as.numeric(dev.list()[length(dev.list())]);
    if(length(v)==0)v=1;
return(v);
}


#' @export
cpal <- function(palid=0,bgid=-2){
        if(class(palid)=="character"){
            a$pparam$resi$bgid <<- -1;
            a$pparam$resi$bg <<- palid;
            cat("\n bg color: ",a$pparam$resi$bg,"\n");
            
        }else{
            palname=c("1: ranbow", "2: heat.colors", "3: terrain.colors", "4: topo.colors", "5: cm.colors");
##            palname.squash=c('6: rainbow2', '7: jet', '8: grayscale', '9: heat', '10: coolheat', '11: blueorange', '12: bluered', '13: darkbluered');
##            palname=c(palname,palname.squash);
            
            
            if((palid==0)&&(bgid==-2)){
            cat("\nPalettes \n", paste(palname,collapse="\n "), "\n");
            }else{
                if(bgid==-2){
                    if(a$pparam$resi$bgid==-1)a$pparam$resi$bgid<<-1;
                }else{
                    a$pparam$resi$bgid<<-bgid;
                }
               
                if(palid==1)a$pparam$pal <<- rev(rainbow(100,end=0.7));
                if(palid==2)a$pparam$pal <<- heat.colors(100);
                if(palid==3)a$pparam$pal <<- terrain.colors(100);
                if(palid==4)a$pparam$pal <<- topo.colors(100);
                if(palid==5)a$pparam$pal <<- cm.colors(100);
                ##if(palid==6)a$pparam$pal <<- rainbow2(100);
                ##if(palid==7)a$pparam$pal <<- jet(100);
                ##if(palid==8)a$pparam$pal <<- grayscale(100);
                ##if(palid==9)a$pparam$pal <<- heat(100);
                ##if(palid==10)a$pparam$pal <<- coolheat(100);
                ##if(palid==11)a$pparam$pal <<- blueorange(100);
                ##if(palid==12)a$pparam$pal <<- bluered(100);
                ##if(palid==13)a$pparam$pal <<- darkbluered(100);
                cat("\n Pallete: ",palname[palid], " bgid:",a$pparam$resi$bgid,"\n");
            }    
        }
}

        
#' @export
pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=600,outeps=FALSE,prefix="./",show=TRUE,outpng=TRUE){
    
    dens=as.character(density);
    geom=as.character(as.integer(geometry));
    dev.set(dev_id);

    tmpf=paste0(a$mydir,"R_pngout_temp");
    
    dev.copy2eps(file=sprintf("%s.eps",tmpf));
   
     if(outeps){
         ##system(paste("cp .R_pngout_temp.eps ",prefix,plotfile,".eps",sep=""));
         system(sprintf("cp %s.eps %s%s.eps",tmpf,prefix,plotfile));
         cat("eps output:",paste(prefix,plotfile,".eps",sep=""),"\n")
     }
    
    if(outpng){
        system(sprintf("convert -density %sx%s -geometry %s  -background white -alpha remove %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        system(sprintf("mv %s.png %s%s.png",tmpf,prefix,plotfile));
        cat(sprintf("png output: %s%s.png \n",prefix,plotfile));
        if(show)system(sprintf("display %s%s.png&",prefix,plotfile));
    
    }
    system(sprintf("rm %s.eps",tmpf));
    
}




#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/ MAIN FUNCTION FOR SIMULATION OF ADAPTIVE EVOLUTION _/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
##' This function simulates adaptive evolution by means of the oligomorphic stochastic model (OSM). The OSM assumes very rare mutations in comparison with the time scale of population dynamics, so that evolutionary dynamics can be described as a trait substitution sequence, engendered by repeated mutant invasions.
##'
##' The population dynamics triggerred by each mutant invasion is calculated with R-package "deSolve".
##' @title simulates adaptive evolution under given fitness function and mutation function.
##' @param phe list: phenotypes of coexisting residents.
##' @param en array: population densities of coexisting residents and environmental variables.
##' @param fitness function: fitness function.
##' @param mutate function: function for mutation.
##' @param pop_dynamics function: function for population dynamics. When this function is not given, the "fitness" function is used for calculation of population dynamics.
##' @return no direct output (the variable "a" contains all simulation data as well as parameters)
##' @author Hiroshi C. Ito
#' @export
simevol <- function(phe=a$phe,en=a$en,## state values
                    fitness=NULL,## functions
                    mutate=mutate0,
                    pop_dynamics=pop_dynamics0,
                    set_parms=NULL,
                    plot_func=plot_func0,
                    output=output0,
                    output_tree=output_tree0,
                    halt_func=NULL,
                    yorick_plot=NULL,
                    tmax=100000000, ## parameters for simulation and output
                    out_interval=10,
                    show_interval=10,
                    file_data="test.dat",
                    file_data_tree="test_tree.dat",
                    file_data_pid="test_pid.dat",
                    continue=FALSE,
                    runid=1,
                    runname="",
                    mydir=".simevol/",
                    fitness_contour=TRUE,## parameters for plotting
                    fitness_contour_phe=NULL,
                    plot_mask=NULL,
                    reset_win=TRUE,
                    show_subwin=TRUE,
                    trait_names=NULL,
                    env_names=NULL,
                    bgid=-1,
                    palid=1,
                    palfunc=NULL,
                    pparam=pparam0,
                    amp_invf=0.1,## parameters for mutant invasion
                    level_invf=0.02,
                    amp_invf_fix=FALSE,
                    m_rate=1.0,
                    m_sd=0.01,
                    n_mutant_init=1e-6,## parameters for population dynamics
                    edge_extinct=1e-8,
                    edge_fit=1e-13,
                    nrad=12,
                    drad=1.0,
                    divrad=2,
                    param_desolve=list(method="radau",hini=1e-4,hmax=1e9,rtol=1e-4,atol=1e-20)                      
                  ){


    ##    .ee.append("adsim",environment())
 
    if(continue==FALSE){ 

        if(!file.exists(mydir))system(sprintf("mkdir %s",mydir));
        file_outcount=paste0(mydir,"outcount.dat");
        file_command=paste0(mydir,"command",runid,".R");
        
        timen=0.0;        
        outcount=1;
        
        nspe=length(phe[[1]]);
        n=en[1:nspe];
        pdim=length(phe);
        edim=length(en)-nspe;
        if(length(trait_names)==0)trait_names=names(phe);
        
        .simevol_func <<- list(fitness=fitness,pop_dynamics=pop_dynamics,mutate=mutate,set_parms=set_parms,plot_func=plot_func,output=output,output_tree=output_tree,halt_func=halt_func,palfunc=palfunc);
        
     
        pparam$plot_mask=plot_mask;
        pparam$trait_names=trait_names;
        pparam$env_names=env_names;
        pparam$fitness_contour=fitness_contour;
        pparam$fitness_contour_phe=fitness_contour_phe;
        pparam$resi$bgid=bgid;

        pparam=c(pparam,list(sub_xlim=vector("list",length=pdim),sub_ylim=vector("list",length=pdim)));

        
        traj=list(phe=phe,n=c(n),t=c(rep(0.0,nspe)),e=c(NULL),te=c(0.0));
        phe=c(phe,list(t=0.0,pid=0));
        tree=list(phe);
        tree_phe=c(phe,list(pid_par=-1));

        
        sparam=list(m_rate=m_rate,
                    m_sd=m_sd,
                invf=c(0.1),
                fit_over=c(NULL),
                flag_halt=FALSE,
                edge_extinct=edge_extinct,
                edge_fit=edge_fit,
                amp_invf=amp_invf,
                level_invf=level_invf,
                amp_invf_fix=amp_invf_fix,
                n_mutant_init=n_mutant_init,
                mydir=mydir,
                file_command=file_command,
                file_outcount=file_outcount,
                outcount=outcount,
                out_interval=out_interval,
                show_interval=show_interval,
                nrad=nrad,
                drad=drad,
                divrad=divrad,
                param_desolve=param_desolve
                );

   

    a<<-list(
        phe=phe,
        en=en,
        n=n,
        e=en[(nspe+1):length(en)],
        timen=timen,
        ninv=0,
        inv=c(NULL),
        edim=edim,
        pdim=pdim,
        traj=traj,
        tree=tree,
        tree_phe=tree_phe,
        winid=c(2,3,0),
        sparam=sparam,
        pparam=pparam,
        file_data=file_data,
        file_data_tree=file_data_tree,
        file_data_pid=file_data_pid,

        runid=runid,
        runname=runname,
        show_subwin=show_subwin
        );
        options(scipen=100);
        write(c(1,0),file=a$file_data_pid,append=FALSE,ncolumns=2);
        options(scipen=0);
        ##cat("\n",file=a$file_data_pid,append=TRUE);


        cpal(palid);
        
        ##comwin();   

        
        if(a$edim>0)a$traj$e<<-c(a$traj$e,en[(nspe+1):length(en)]);
        
        if(file.exists(file_data))file.remove(file_data);
        
        .simevol_func$output(a$timen,phe,n);
        cat(a$sparam$outcount,file=a$sparam$file_outcount);
        a$sparam$outcount<<-a$sparam$outcount+1;
        
        if(class(yorick_plot)=="function"){
            yorick_plot();
            system("xterm -fn 7x14 -bg navy -fg white -e 'rlwrap -c ./yorick_idl_follow.sh'&");
        }
    
        if(reset_win || (length(dev.list())<2)){
            graphics.off();
            plot_init();
        }

        res=simpop(phe,en,.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad);
        

        phe=res$phe;
        nspe=length(phe[[1]]);
        en=res$en;
        n=en[1:nspe];


        if(file.exists(file_command))file.remove(file_command);
    }
    else{
        a$sparam$flag_halt<<-FALSE;
    }

    
   
    for(t in 1:tmax){
        if(class(.simevol_func$halt_func)=="function").simevol_func$halt_func();
        if(a$sparam$flag_halt==T)break;

        if(file.exists(file_command)){
            ##source(file_command);
            sys.source(file_command,envir=.GlobalEnv);
            file.remove(file_command);
        }
        
        
        inv=invader(phe,en,mutate=.simevol_func$mutate,a$timen,amp_invf=a$sparam$amp_invf);
        
        a$timen<<-inv$timen;
        a$sparam$invf<<-c(a$sparam$invf,inv$f);
        a$inv <<- inv;
        a$ninv <<- a$ninv+1;
        
        state_new=simpop_invade(phe,en,inv,.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad);

        
        par_alive= which(state_new$phe$pid==inv$pid_par);
        a$tree_phe <<- add_phenotype(a$tree_phe,c(inv$phe,list(pid_par=inv$pid_par)));

        options(scipen=100);
        write(c(inv$phe$pid+1,inv$pid_par+1),file=a$file_data_pid,append=TRUE,ncolumns=2); ## to be improved so that pids in simevol and file_data_pid are the same.
        options(scipen=0);
                
        
        if(length(par_alive)>0){
            ##print(inv$phe$pid);
            ##print(state_new$phe$pid);            
            ##print(par_alive);
            
            phe_par=geti(state_new$phe,par_alive);
            ##print(phe_par);
            a$tree <<- c(a$tree,list(add_phenotype(phe_par,inv$phe)));
                
        }
        else{
            for(i in 1:length(a$tree)){
                buf=a$tree[[i]];
                ##cat("pars:",buf$pid);
                ##cat("inv:",inv$phe$pid,"inv_par:",inv$pid_par,"\n");
                if(sum(buf$pid[length(buf$pid)]==inv$pid_par)>0){
                    ##cat("connect\n");
                    a$tree[[i]]<<-add_phenotype(a$tree[[i]],inv$phe);
                }
            }
        }
        
        for(i in 1:length(a$tree)){
            buf=a$tree[[i]];
            buf=geti(buf,length(buf$pid));
            if(buf$pid>0){
                if(sum(state_new$phe$pid==buf$pid)==0){
                    ##print(buf$pid);
                    ##   a$tree_phe$te[(buf$pid+1)]<<-a$timen;
                    buf$t=a$timen;
                    buf$pid=-1; ## -1 means extinction
                    a$tree[[i]]<<-add_phenotype(a$tree[[i]],buf);
                }
                
            }
        }
        
   

        
        phe=state_new$phe;
        nspe=length(phe[[1]]);
        en=state_new$en;
        n=en[1:nspe];
            
        a$phe<<-phe;
        a$en<<-en;
        a$n<<-n;
        a$e<<- en[(nspe+1):length(en)];


        level_invf=a$sparam$level_invf;
        if(amp_invf_fix==FALSE){
            a$sparam$amp_invf<<-max(level_invf,level_invf/mean(a$sparam$invf[max((t-100),1):length(a$sparam$invf)]));
        }

        
        .simevol_func$output(a$timen,phe,n);
                
        cat(a$sparam$outcount,file=a$sparam$file_outcount);
        a$sparam$outcount<<-a$sparam$outcount+1;

        
        if((a$sparam$flag_halt==T)||(a$ninv%%a$sparam$out_interval==0)){
            a$traj$phe<<-add_phenotype(a$traj$phe,phe);
            .simevol_func$output_tree(a$file_data_tree);
        
            
            if(a$edim>0)a$traj$e<<-rbind(a$traj$e,en[(nspe+1):length(en)]);
            a$traj$n<<-c(a$traj$n,n);
            a$traj$t<<-c(a$traj$t,rep(a$timen,nspe));
            a$traj$te<<-c(a$traj$te,a$timen);
            
        }
    
        if((a$sparam$flag_halt==T)||(a$ninv%%a$sparam$show_interval==0)){
            runname1="";
            if(runname!="")runname1=paste0("\"",runname,"\"");
            cat("runid:",runid,runname1,"time:",a$timen, " residents:",nspe, " invasion:", a$ninv,"amp:",a$sparam$amp_invf,"fit_over:",length(a$sparam$fit_over),"irad:",state_new$irad,"\n");
            
            .simevol_func$plot_func();    
            
        }
##        if(a$sparam$flag_halt==T)break;
    }
    
}


