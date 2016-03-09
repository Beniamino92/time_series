setwd("Module_8_Time_Series/data/")
ts_data <- lapply(list.files(), function(x) read.csv(x))

y <- ts_data[[6]]
y[,2] <- as.POSIXct(y[,2])
y[,3] <- (strftime(y[,2], format = "%H") < 12)
ggplot(y, aes(DateTime, Value, V3)) +
  geom_rect(aes(xmin=DateTime-1800, xmax=DateTime+1800, 
                ymin=min(Value), ymax=max(Value), 
                fill=V3), show.legend = F) + 
  scale_fill_manual(values = c("#7EC0EE33", "#00007c33"))+
  geom_line()+geom_point()+
  scale_x_datetime(date_breaks = "day", expand = c(0,0) ) + xlab("") + ylab("Rest Activity")+
  scale_y_continuous(expand = c(0,0))+theme_bw()



