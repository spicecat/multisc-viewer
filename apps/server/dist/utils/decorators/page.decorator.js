"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.PAGE_METADATA = void 0;
exports.Page = Page;
exports.PAGE_METADATA = Symbol("PAGE");
function Page() {
    return (target, key, descriptor) => {
        Reflect.defineMetadata(exports.PAGE_METADATA, true, descriptor.value);
        return descriptor;
    };
}
//# sourceMappingURL=page.decorator.js.map