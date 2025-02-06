import type { ArgumentsHost } from "@nestjs/common";
import type { Response } from "express";
import { Catch, HttpException } from "@nestjs/common";
import { BaseExceptionFilter } from "@nestjs/core";

export class ErrorPage {
  public constructor(public readonly props: any) {}
}

@Catch()
export class ErrorPageFilter extends BaseExceptionFilter {
  public catch(exception: unknown, host: ArgumentsHost): void {
    const ctx = host.switchToHttp();
    const response = ctx.getResponse<Response>();

    if (exception instanceof HttpException) {
      const err = exception.getResponse();

      if (err instanceof ErrorPage) {
        response.status(exception.getStatus()).render("error", err.props);
        return;
      }
    }

    super.catch(exception, host);
  }
}
