#!/usr/bin/Rscripts
data.df <-
    read.csv("Power_Data_30_min_kWh.csv", stringsAsFactors=FALSE);

colnames(data.df) <- c("Date","Time","Usage");

data.df$Date <- as.POSIXct(data.df$Date, format="%d/%m/%Y");
data.df$Time <- as.POSIXct(data.df$Time, format="%I:%M:%S %p");

## Graph of power usage by half-hour (aggregated over all days), with median
png("PowerUsageDayAggregate.png", width=1280, height=640, pointsize=24);
par(xaxs="i", mar=c(4,4,0.5,0.5), cex.axis=0.8);
smoothScatter(y=log10(data.df$Usage),
              x= as.numeric(data.df$Time),
     pch=21, bg="#00000010",
     col="#00000010", xaxt="n", yaxt="n", ylab = "Usage (kWh)",
     xlab = "Time", nbin=c(48*4,128), bandwidth=c(1200,.04), nrpoints=0);
axis(1, at=sort(unique(data.df$Time)), cex.axis=0.71, las=2,
     labels=strftime(sort(unique(data.df$Time)), format="%H:%M %p"));
axis(2, at=log10(1:9 * rep(10^(-1:1),each=9)),
     labels=(1:9 * rep(10^(-1:1),each=9)), las=1);
day.aggregate <- tapply(data.df$Usage, data.df$Time, median);
points(spline(x=as.POSIXct(names(day.aggregate)),y=log10(day.aggregate)),
       type = "l");
points(x=as.POSIXct(names(day.aggregate)),y=log10(day.aggregate), lwd=2);
invisible(dev.off());

## Graph of daily power usage, with median daily usage per month
png("PowerUsageByMonth.png", width=1280, height=640, pointsize=24);
par(mar=c(4,4,0.5,0.5), cex.axis=0.8);
plot(x=as.POSIXct(names(day.sum)), y=day.sum, xaxt="n", xlab="",
     ylab = "Daily Consumption (kWh)", col="#00000020", pch=21,
     bg="#00000020", cex=0.5);
axis(1, at=as.POSIXct(sprintf("201%d-%02d-01",rep(5:7,each=11), 2:12)),
     label=strftime(sprintf("201%d-%02d-01",rep(5:7,each=11), 2:12),
                    format="%b"), las=2);
axis(1, at=as.POSIXct(sprintf("201%d-%02d-01",rep(5:7,each=1), 1)),
     label=strftime(sprintf("201%d-%02d-01",rep(5:7,each=1), 1),
                    format="%Y-%b"), las=2);
day.sum <- tapply(data.df$Usage, data.df$Date, sum);
month.med <- tapply(day.sum, strftime(names(day.sum), format="%Y-%b"), median);
names(month.med) <- sub("$","-15",names(month.med));
month.med <- month.med[order(as.POSIXct(names(month.med), format="%Y-%b-%d"))];
points(x=as.POSIXct(names(month.med), format="%Y-%b-%d"), y=month.med,
       pch=21, bg="black");
points(spline(x=as.POSIXct(names(month.med), format="%Y-%b-%d"), y=month.med,
              n=length(month.med) * 10),
       type="l");
invisible(dev.off());
