import type { NestExpressApplication } from "@nestjs/platform-express";
import { ConfigService } from "@nestjs/config";

import { ValidationPipe } from "@nestjs/common";
import { HttpAdapterHost, NestFactory } from "@nestjs/core";
// import cookieParser from "cookie-parser";
import { AppModule } from "./app.module";
// import { svelte } from "./client/template-engine.js";
// import { ErrorPageFilter } from "./utils/filters/error-page.filter";
// import { RedirectFilter } from "./utils/filters/redirect.filter";
// import { RoutingInterceptor } from "./utils/interceptors/routing.interceptor";

async function bootstrap() {
  const app = await NestFactory.create<NestExpressApplication>(AppModule);
  const configService = app.get(ConfigService);
  // app.engine("svelte", svelte);
  // app.setViewEngine("svelte");
  // app.setBaseViewsDir("src/client/routes");

  // app
  // .use(cookieParser())
  // .useGlobalPipes(
  //   new ValidationPipe({
  //     transform: true,
  //     transformOptions: { enableImplicitConversion: true },
  //   })
  // )
  // .useGlobalFilters(
  //   new ErrorPageFilter(app.get(HttpAdapterHost).httpAdapter),
  //   new RedirectFilter()
  // )
  // .useGlobalInterceptors(new RoutingInterceptor());

  await app.listen(configService.get<number>('app.port', 5000));
}

bootstrap();
