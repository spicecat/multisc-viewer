"use strict";
var __decorate = (this && this.__decorate) || function (decorators, target, key, desc) {
    var c = arguments.length, r = c < 3 ? target : desc === null ? desc = Object.getOwnPropertyDescriptor(target, key) : desc, d;
    if (typeof Reflect === "object" && typeof Reflect.decorate === "function") r = Reflect.decorate(decorators, target, key, desc);
    else for (var i = decorators.length - 1; i >= 0; i--) if (d = decorators[i]) r = (c < 3 ? d(r) : c > 3 ? d(target, key, r) : d(target, key)) || r;
    return c > 3 && r && Object.defineProperty(target, key, r), r;
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.RoutingInterceptor = void 0;
const common_1 = require("@nestjs/common");
const constants_1 = require("@nestjs/common/constants");
const node_fs_1 = require("node:fs");
const rxjs_1 = require("rxjs");
const operators_1 = require("rxjs/operators");
const page_decorator_1 = require("../decorators/page.decorator");
function join(route, part) {
    return route === "" ? part : `${route}/${part}`;
}
function resolveRoute(parts, route, dir, params) {
    if (parts.length === 0) {
        return route;
    }
    const available = (0, node_fs_1.readdirSync)(dir);
    let pattern;
    if (parts.length === 1) {
        const [part] = parts;
        if (available.includes(`${part}.svelte`)) {
            return join(route, part);
        }
        else if (part === "" && available.includes("index.svelte")) {
            return join(route, "index");
        }
        else if (available.includes(part) &&
            (0, node_fs_1.statSync)(`${dir}/${part}`).isDirectory() &&
            (0, node_fs_1.readdirSync)(`${dir}/${part}`).includes("index.svelte")) {
            return join(route, `${part}/index`);
        }
        else if ((pattern = available.find((route) => {
            const match = /\[(\w+)\]\.svelte/.exec(route);
            return match !== null && match[1] in params;
        })) !== undefined) {
            return join(route, pattern.replace(".svelte", ""));
        }
        else {
            return null;
        }
    }
    else {
        const [part, ...rest] = parts;
        if (available.includes(part) && (0, node_fs_1.statSync)(`${dir}/${part}`).isDirectory()) {
            return resolveRoute(rest, join(route, part), `${dir}/${part}`, params);
        }
        else if ((pattern = available.find((route) => {
            const match = /\[(\w+)\]/.exec(route);
            return match !== null && match[1] in params;
        })) !== undefined &&
            (0, node_fs_1.statSync)(`${dir}/${pattern}`).isDirectory()) {
            return resolveRoute(rest, join(route, pattern), `${dir}/${pattern}`, params);
        }
        else {
            return null;
        }
    }
}
let RoutingInterceptor = class RoutingInterceptor {
    intercept(context, next) {
        const req = context.switchToHttp().getRequest();
        const res = context.switchToHttp().getResponse();
        // this case is mostly for backwards compatibility
        if (Reflect.hasMetadata(constants_1.RENDER_METADATA, context.getHandler())) {
            const path = Reflect.getMetadata(constants_1.RENDER_METADATA, context.getHandler());
            const parts = path
                .split("/")
                .map((part, i, arr) => i === arr.length - 1 ? part.replace(".svelte", "") : part);
            const route = parts.slice(parts.indexOf("routes") + 1).join("/");
            return next.handle().pipe((0, operators_1.map)((val) => {
                const props = val ?? {};
                const user = props.user ?? req.user;
                const defaultMeta = {
                    route,
                    path: req.path,
                    user,
                    params: req.params,
                    query: req.query,
                };
                const __meta = props.__meta
                    ? { ...defaultMeta, ...props.__meta }
                    : defaultMeta;
                if ("user" in props)
                    delete props.user;
                if ("__meta" in props)
                    delete props.__meta;
                return { props, __meta };
            }));
        }
        else if (Reflect.hasMetadata(page_decorator_1.PAGE_METADATA, context.getHandler())) {
            const parts = req.path.split("/");
            const route = resolveRoute(parts.slice(1), "", "src/client/routes", req.params);
            if (route !== null) {
                return next.handle().pipe((0, operators_1.map)((val) => {
                    const props = val ?? {};
                    const user = props.user ?? req.user;
                    const defaultMeta = {
                        route,
                        path: req.path,
                        user,
                        params: req.params,
                        query: req.query,
                    };
                    const __meta = props.__meta
                        ? { ...defaultMeta, ...props.__meta }
                        : defaultMeta;
                    if ("user" in props)
                        delete props.user;
                    if ("__meta" in props)
                        delete props.__meta;
                    return (0, rxjs_1.from)(new Promise((resolve, reject) => res.render(route, { props, __meta }, (err, html) => err !== null ? reject(err) : resolve(html))));
                }), (0, operators_1.mergeAll)());
            }
            else {
                return (0, rxjs_1.throwError)(() => new common_1.NotFoundException()); // TODO: consider adding message
            }
        }
        else {
            return next.handle();
        }
    }
};
exports.RoutingInterceptor = RoutingInterceptor;
exports.RoutingInterceptor = RoutingInterceptor = __decorate([
    (0, common_1.Injectable)()
], RoutingInterceptor);
//# sourceMappingURL=routing.interceptor.js.map