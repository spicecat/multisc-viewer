import type { ArgumentsHost, ExceptionFilter } from "@nestjs/common";
import type { Response } from "express";
import { Catch, HttpStatus } from "@nestjs/common";

export class Redirect {
  public constructor(
    public readonly location: string,
    public readonly status: number = HttpStatus.SEE_OTHER,
  ) {}
}

@Catch(Redirect)
export class RedirectFilter implements ExceptionFilter<Redirect> {
  public catch(exception: Redirect, host: ArgumentsHost) {
    const ctx = host.switchToHttp();
    const response = ctx.getResponse<Response>();

    response
      .status(exception.status)
      .setHeader("Location", exception.location)
      .end();
  }
}
