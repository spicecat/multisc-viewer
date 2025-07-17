import { Module } from "@nestjs/common";
import { ConfigModule } from "@nestjs/config";
// import { ServeStaticModule } from "@nestjs/serve-static";
import { AppController } from "./app.controller";
import configuration from "./config/configuration";
import { DaemonService } from "./daemon.service";
import { DataService } from "./data.service";
import { PlotService } from "./plot.service";

@Module({
  imports: [
    ConfigModule.forRoot({
      load: [configuration],
    }),
    // ServeStaticModule.forRoot({
    //   rootPath: "dist/client/assets",
    //   serveRoot: "/__app",
    // }),
    // ServeStaticModule.forRoot({
    //   rootPath: "src/client/public",
    //   serveRoot: "/",
    // }),
  ],
  controllers: [AppController],
  providers: [DataService, PlotService, DaemonService],
})
export class AppModule {}
