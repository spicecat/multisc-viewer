import { Module } from "@nestjs/common";
import { ConfigModule } from "@nestjs/config";
// import { ServeStaticModule } from "@nestjs/serve-static";
import configuration from "./config/configuration";
import { AppController } from "./app.controller";
import { AppService } from "./app.service";
import { DaemonService } from "./daemon.service";

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
  providers: [AppService, DaemonService],
})
export class AppModule {}
